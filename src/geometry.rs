use std::{
    error::Error,
    path::Path,
    convert::TryInto,
};

use geo::{
    algorithm::area::Area,
    BoundingRect,
    Geometry,
    MultiPolygon,
    Coord,
    LineString,
    MapCoords,
    Polygon,
};
use geos::{
    Geometry as GeosGeometry,
    Geom,
};

use rstar::{AABB, RTree};
use rayon::prelude::*;
use ndarray::{s, Array, Array2, Zip};
use proj::Proj;

use crate::utils::{RvnGridWeights, AREA_THRESHOLD};
use crate::io::read_lat_lon;
use crate::cli::Cli;
use crate::rtree::{build_grid_rtree, GridCell};


// Create a meshgrid from 1D arrays of latitude and longitude
pub fn meshgrid(lat: &Array<f32, ndarray::Ix1>, lon: &Array<f32, ndarray::Ix1>) -> (Array<f32, ndarray::Ix2>, Array<f32, ndarray::Ix2>) {
    // Create 2D longitude by repeating lon along the rows
    let lon_2d = Array::from_shape_fn((lat.len(), lon.len()), |(_, j)| lon[j]);

    // Create 2D latitude by repeating lat along the columns
    let lat_2d = Array::from_shape_fn((lat.len(), lon.len()), |(i, _)| lat[i]);

    (lat_2d, lon_2d)
}

// Creates the coordinates of the grid cells using the centers of the grid cells
pub fn create_gridcells_from_centers(lat: &Array2<f32>, lon: &Array2<f32>) -> (Array2<f32>, Array2<f32>) {
    let nlat = lat.nrows();
    let nlon = lon.ncols();

    // Create arrays for the edges
    let mut lonh = Array2::<f32>::zeros((nlat + 1, nlon + 1));
    let mut lath = Array2::<f32>::zeros((nlat + 1, nlon + 1));

    // Temporary arrays for dlat and dlon
    let mut dlat = Array2::<f32>::zeros((nlat, nlon));
    let mut dlon = Array2::<f32>::zeros((nlat, nlon));

    // Calculate dlat and dlon
    for ii in 0..(nlat - 1) {
        for jj in 0..(nlon - 1) {
            dlat[[ii, jj]] = (lat[[ii + 1, jj + 1]] - lat[[ii, jj]]) / 2.0;
            dlon[[ii, jj]] = (lon[[ii + 1, jj + 1]] - lon[[ii, jj]]) / 2.0;
        }
        dlat[[ii, nlon - 1]] = (lat[[ii + 1, nlon - 1]] - lat[[ii, nlon - 2]]) / 2.0;
        dlon[[ii, nlon - 1]] = (lon[[ii + 1, nlon - 1]] - lon[[ii, nlon - 2]]) / 2.0;
    }

    // Duplicate the last row of dlat and dlon
    let last_row_dlat = dlat.row(nlat - 2).to_owned();
    let last_row_dlon = dlon.row(nlat - 2).to_owned();
    dlat.slice_mut(s![nlat - 1, ..]).assign(&last_row_dlat);
    dlon.slice_mut(s![nlat - 1, ..]).assign(&last_row_dlon);

    // Compute lath and lonh by adjusting the lat and lon values by their respective deltas
    Zip::from(&mut lonh.slice_mut(s![..nlat, ..nlon]))
        .and(&*lon)
        .and(&dlon)
        .for_each(|lonh_elem, &lon_val, &dlon_val| {
            *lonh_elem = lon_val - dlon_val;
        });

    Zip::from(&mut lath.slice_mut(s![..nlat, ..nlon]))
        .and(&*lat)
        .and(&dlat)
        .for_each(|lath_elem, &lat_val, &dlat_val| {
            *lath_elem = lat_val - dlat_val;
        });

    // Extend the grid at the edges
    for jj in 0..nlon {
        lonh[[nlat, jj]] = lonh[[nlat - 1, jj]] + (lonh[[nlat - 1, jj]] - lonh[[nlat - 2, jj]]);
        lath[[nlat, jj]] = lath[[nlat - 1, jj]] + (lath[[nlat - 1, jj]] - lath[[nlat - 2, jj]]);
    }

    for ii in 0..nlat {
        lonh[[ii, nlon]] = lonh[[ii, nlon - 1]] + (lonh[[ii, nlon - 1]] - lonh[[ii, nlon - 2]]);
        lath[[ii, nlon]] = lath[[ii, nlon - 1]] + (lath[[ii, nlon - 1]] - lath[[ii, nlon - 2]]);
    }

    lonh[[nlat, nlon]] = lonh[[nlat - 1, nlon - 1]] + (lonh[[nlat - 1, nlon - 1]] - lonh[[nlat - 2, nlon - 2]]);
    lath[[nlat, nlon]] = lath[[nlat - 1, nlon - 1]] + (lath[[nlat - 1, nlon - 1]] - lath[[nlat - 2, nlon - 2]]);

    (lath, lonh)
}

// Convert a list of coordinates to a Polygon
pub fn shape_to_geometry(
    shape_coord: &[[f64; 2]],
    proj: &Proj,
) -> Result<Polygon<f64>, Box<dyn Error>> {
    // Convert the list of coordinates to a Vec of Coordinates
    let mut coords: Vec<Coord<f64>> = shape_coord
        .iter()
        .map(|&point| Coord {
            x: point[1], // Swap lat/lon to x/y
            y: point[0],
        })
        .collect();

    // Close the ring by adding the first point at the end if necessary
    if coords.first() != coords.last() {
        coords.push(coords[0]);
    }

    // Before projection, adjust latitudes at the poles
    for coord in &mut coords {
        if coord.y > 90.0 {
            coord.y = 89.9999;
        } else if coord.y < -90.0 {
            coord.y = -89.9999;
        }
    }

    let mut exterior = LineString::from(coords);
    exterior = exterior
        .try_map_coords(|coord| proj.convert(coord)).expect("Failed to transform coordinates");

    let polygon = Polygon::new(exterior, vec![]);

    Ok(polygon)
}

/// Process lat/lon and create grid cells
pub fn process_lat_lon(args: &Cli) -> Result<(Array2<f32>, Array2<f32>, usize, usize), Box<dyn Error>> {
    let file = netcdf::open(Path::new(&args.nc)).expect("Failed to open file.");
    let (lat_2d, lon_2d) = read_lat_lon(&file, &args).expect("Failed to read lat/lon");

    let (lath, lonh) = if args.grd_bnds {
        (lat_2d, lon_2d)
    } else {
        create_gridcells_from_centers(&lat_2d, &lon_2d)
    };
    
    // -1 here because we are going to have 1 polygon between 2 latitudes/longitudes
    let nlat = lath.nrows() - 1;
    let nlon = lonh.ncols() - 1;

    Ok((lath, lonh, nlat, nlon))
}

/// Create grid cells and store them as geometries 
/// TODO: Don't use Vec of Vecs and use a single Vec instead.
pub fn create_grid_cells(
    nlat: usize, 
    nlon: usize, 
    lath: &Array2<f32>, 
    lonh: &Array2<f32>,
    target_epsg: &str,
) -> Result<Vec<Vec<Polygon<f64>>>, Box<dyn Error>> {
    let mut grid_cell_geom: Vec<Vec<Polygon<f64>>> = Vec::with_capacity(nlat);
    let proj = Proj::new_known_crs("EPSG:4326", target_epsg, None).expect("Failed to create projection");

    for ilat in 0..nlat {
        let mut row: Vec<Polygon<f64>> = Vec::with_capacity(nlon);
        for ilon in 0..nlon {
            let gridcell_edges = [
                [lath[[ilat, ilon]] as f64, lonh[[ilat, ilon]] as f64],
                [lath[[ilat + 1, ilon]] as f64, lonh[[ilat + 1, ilon]] as f64],
                [lath[[ilat + 1, ilon + 1]] as f64, lonh[[ilat + 1, ilon + 1]] as f64],
                [lath[[ilat, ilon + 1]] as f64, lonh[[ilat, ilon + 1]] as f64],
            ];
            let polygon = shape_to_geometry(&gridcell_edges, &proj).expect("Failed to create geometry");
            row.push(polygon);
        }
        grid_cell_geom.push(row);
    }
    Ok(grid_cell_geom)
}

/// Check if the bounding boxes intersect by comparing their coordinates
fn bounding_boxes_intersect(bbox1: &geo_types::Rect<f64>, bbox2: &geo_types::Rect<f64>) -> bool {
    !(bbox1.max().x < bbox2.min().x || bbox1.min().x > bbox2.max().x || 
      bbox1.max().y < bbox2.min().y || bbox1.min().y > bbox2.max().y)
}

/// Calculate the intersection between the grid cell and the shape geometry
fn calculate_intersection(grid_cell: &Polygon<f64>, shape_geometry: &Geometry) -> Option<Geometry> {
    // Check the bounding boxes of the grid cell and shape geometry
    let grid_cell_bbox = grid_cell.bounding_rect()?;
    let shape_bbox = shape_geometry.bounding_rect()?;
    if !bounding_boxes_intersect(&grid_cell_bbox, &shape_bbox) {
        return None;
    }

    // Convert grid cell to GEOS Geometry and calculate the intersection
    let grid_cell_geos: GeosGeometry = grid_cell.try_into().expect("Failed to convert Polygon to GEOS Geometry");
    match shape_geometry {
        Geometry::Polygon(ref shape_poly) => {
            // Convert shape polygon to GEOS Geometry
            let shape_poly_geos: GeosGeometry = shape_poly.try_into().ok()?;
            let intersection = grid_cell_geos.intersection(&shape_poly_geos).ok()?;

            if !intersection.is_empty().ok()? {
                let result: Geometry<f64> = intersection.try_into().ok()?;
                Some(result)
            } else {
                None
            }
        },
        Geometry::MultiPolygon(ref shape_mpoly) => {
            // Convert MultiPolygon to GEOS Geometry
            let shape_mpoly_geos: GeosGeometry = shape_mpoly.try_into().ok()?;
            let intersection = grid_cell_geos.intersection(&shape_mpoly_geos).ok()?;

            if !intersection.is_empty().ok()? {
                let result: Geometry<f64> = intersection.try_into().ok()?;
                Some(result)
            } else {
                None
            }
        },
        _ => None,
    }
}

/// Calculate the area of a geometry
fn calculate_area(geometry: &Geometry) -> f64 {
    match geometry {
        Geometry::Polygon(ref poly) => poly.unsigned_area(),
        Geometry::MultiPolygon(ref mpoly) => mpoly.unsigned_area(),
        _ => 0.0, // Since it can be a LineString.
    }
}

/// Used to process a single shape and calculate the intersection with grid cells
fn process_shape(
    shape_geometry: &Geometry,
    rtree: &RTree<GridCell>,
    nlon: usize,
) -> Option<(Vec<(usize, f64)>, f64, f64)> {
    // Calculate the area of the shape
    let shape_area = match shape_geometry {
        Geometry::Polygon(ref poly) => poly.unsigned_area(),
        Geometry::MultiPolygon(ref mpoly) => mpoly.unsigned_area(),
        _ => return None, // Skip if not a polygon or multipolygon
    };

    // Get the bounding rectangle of the shape
    let shape_envelope = shape_geometry.bounding_rect()?;

    let mut area_all = 0.0;
    let mut cellid_fraction = Vec::new();

    // Query the R-tree for grid cells intersecting the shape's bounding rectangle
    let candidate_cells = rtree.locate_in_envelope_intersecting(&AABB::from_corners(
        [shape_envelope.min().x, shape_envelope.min().y],
        [shape_envelope.max().x, shape_envelope.max().y],
    ));

    // Process intersections
    for grid_cell in candidate_cells {
        let intersection = calculate_intersection(&grid_cell.polygon, shape_geometry);
        if let Some(intersection_geom) = intersection {
            let area_intersect = calculate_area(&intersection_geom);
            if area_intersect > 0.0 {
                area_all += area_intersect;
                let cell_id = grid_cell.ilat * nlon + grid_cell.ilon;
                let fraction = area_intersect / shape_area;
                cellid_fraction.push((cell_id, fraction));
            }
        }
    }

    let error = (shape_area - area_all) / shape_area;
    Some((cellid_fraction, area_all, error))
}

fn process_single_shape_netcdf(
    shape_geometry: &Geometry,
    rtree: &RTree<GridCell>,
    nlon: usize,
    area_error_threshold: f64,
    index: usize,
) -> Vec<(f64, usize, usize)> {
    let mut netcdf_data_chunk = Vec::new();

    if let Some((cellid_fraction, _area_all, error)) = process_shape(shape_geometry, rtree, nlon) {
        if error.abs() > area_error_threshold {
            for (cell_id, fraction) in cellid_fraction {
                netcdf_data_chunk.push((fraction, cell_id + 1, index + 1));
            }
        } else {
            let correction_factor = 1.0 / (1.0 - error);
            for (cell_id, fraction) in cellid_fraction {
                let corrected = fraction * correction_factor;
                netcdf_data_chunk.push((corrected, cell_id + 1, index + 1));
            }
        }
    }

    netcdf_data_chunk
}

fn process_single_shape_txt(
    basin_id: &String,
    shape_geometry: &Geometry,
    rtree: &RTree<GridCell>,
    nlon: usize,
    area_error_threshold: f64,
) -> Vec<(String, usize, f64)> {
    let mut txt_data_chunk = Vec::new();

    if let Some((cellid_fraction, _area_all, error)) = process_shape(shape_geometry, rtree, nlon) {
        if error.abs() > area_error_threshold {
            for (cell_id, fraction) in cellid_fraction {
                txt_data_chunk.push((basin_id.clone(), cell_id, fraction));
            }
        } else {
            let correction_factor = 1.0 / (1.0 - error);
            for (cell_id, fraction) in cellid_fraction {
                let corrected = fraction * correction_factor;
                txt_data_chunk.push((basin_id.clone(), cell_id, corrected));
            }
        }
    }

    txt_data_chunk
}

/// Process the shape intersections with grid cells and collect results for netcdf
pub fn process_shape_intersections_netcdf(
    nlon: usize, 
    grid_cell_geom: Vec<Vec<Polygon<f64>>>, 
    shapes: Vec<(String, Geometry)>, 
) -> Result<Vec<(f64, usize, usize)>, Box<dyn Error>> {
    let mut netcdf_data = Vec::new();

    // 5 percent area error threshold. Is 5 percent enough?
    let area_error_threshold = AREA_THRESHOLD;

    // Build the R-tree from grid cells
    let rtree = build_grid_rtree(grid_cell_geom);

    // Iterate over the basins and calculate the intersection with grid cells
    for (index, (_, shape_geometry)) in shapes.into_iter().enumerate() {
        let netcdf_chunk = process_single_shape_netcdf(
            &shape_geometry, &rtree, nlon, area_error_threshold, index,
        );
        netcdf_data.extend(netcdf_chunk);
    }
    Ok(netcdf_data)
}

/// Process the shape intersections with grid cells and collect results for netcdf
pub fn process_shape_intersections_txt(
    nlat: usize, 
    nlon: usize, 
    grid_cell_geom: Vec<Vec<Polygon<f64>>>, 
    shapes: Vec<(String, Geometry)>, 
) -> Result<RvnGridWeights, Box<dyn Error>> {
    let nsubbasins = shapes.len() as i32;
    let ncells = (nlat * nlon) as i32;
    let mut txt_data = Vec::new();

    // 5 percent area error threshold. Is 5 percent enough?
    let area_error_threshold = AREA_THRESHOLD;

    // Build the R-tree from grid cells
    let rtree = build_grid_rtree(grid_cell_geom);

    // Iterate over the basins and calculate the intersection with grid cells
    for (basin_id, shape_geometry) in shapes.into_iter() {
        let txt_chunk = process_single_shape_txt(
            &basin_id, &shape_geometry, &rtree, nlon, area_error_threshold,
        );
        txt_data.extend(txt_chunk);
    }

    let rvn_data = RvnGridWeights {
        txt_data,
        nsubbasins,
        ncells,
    };

    Ok(rvn_data)
}

/// Parallel process the shape intersections with grid cells and collect results
pub fn parallel_process_shape_intersections_netcdf(
    nlon: usize,
    grid_cell_geom: Vec<Vec<Polygon<f64>>>,
    shapes: Vec<(String, Geometry)>,
) -> Result<Vec<(f64, usize, usize)>, Box<dyn Error>> {

    // 5 percent area error threshold
    let area_error_threshold = AREA_THRESHOLD;

    // Build the R-tree from grid cells
    let rtree = build_grid_rtree(grid_cell_geom);

    // Use Rayon to parallelize the loop
    let results: Vec<_> = shapes
        .par_iter()
        .enumerate()
        .map(|(index, (_, shape_geometry))| {
            process_single_shape_netcdf(
                shape_geometry, &rtree, nlon, area_error_threshold, index,
            )
        })
        .collect();

    // Combine the results from all threads
    let mut netcdf_data = Vec::new();

    for netcdf_chunk in results {
        netcdf_data.extend(netcdf_chunk);
    }

    Ok(netcdf_data)
}

/// Parallel process the shape intersections with grid cells and collect results
pub fn parallel_process_shape_intersections_txt(
    nlat: usize,
    nlon: usize,
    grid_cell_geom: Vec<Vec<Polygon<f64>>>,
    shapes: Vec<(String, Geometry)>,
) -> Result< RvnGridWeights, Box<dyn Error>> {
    let nsubbasins = shapes.len() as i32;
    let ncells = (nlat * nlon) as i32;

    // 5 percent area error threshold
    let area_error_threshold = AREA_THRESHOLD;

    // Build the R-tree from grid cells
    let rtree = build_grid_rtree(grid_cell_geom);

    // Use Rayon to parallelize the loop
    let results: Vec<_> = shapes
        .par_iter()
        .enumerate()
        .map(|(_, (basin_id, shape_geometry))| {
            process_single_shape_txt(
                basin_id, shape_geometry, &rtree, nlon, area_error_threshold,
            )
        })
        .collect();

    // Combine the results from all threads
    let mut txt_data = Vec::new();

    for txt_chunk in results {
        txt_data.extend(txt_chunk);
    }

    let rvn_data = RvnGridWeights {
        txt_data,
        nsubbasins,
        ncells,
    };

    Ok(rvn_data)
}


#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{array, Array2};
    use geo::{
        CoordsIter,
        Coord,
        Geometry, 
        polygon, 
        Polygon, 
        algorithm::contains::Contains,
    }; 
    use crate::cli::Cli;
    use std::path::Path;
    use clap::Parser;
    use rstar::RTree;

    #[test]
    fn test_meshgrid() {
        let lat = array![10.0, 20.0, 30.0];
        let lon = array![40.0, 50.0];
        let (lat_2d, lon_2d) = meshgrid(&lat, &lon);
        let expected_lat = array![
            [10.0, 10.0],
            [20.0, 20.0],
            [30.0, 30.0]
        ];
        let expected_lon = array![
            [40.0, 50.0],
            [40.0, 50.0],
            [40.0, 50.0]
        ];
        assert_eq!(lat_2d, expected_lat);
        assert_eq!(lon_2d, expected_lon);
    }

    #[test]
    fn test_create_gridcells_from_centers_square_grid() {
        let lat = array![
            [10.0, 20.0],
            [30.0, 40.0]
        ];
        let lon = array![
            [100.0, 110.0],
            [120.0, 130.0]
        ];
        let (lath, lonh) = create_gridcells_from_centers(&lat, &lon);
        let expected_lath = array![
            [-5.0, 5.0, 15.0],
            [15.0, 25.0, 35.0],
            [35.0, 45.0, 55.0]
        ];
        let expected_lonh = array![
            [85.0, 95.0, 105.0],
            [105.0, 115.0, 125.0],
            [125.0, 135.0, 145.0]
        ];
        assert_eq!(lath, expected_lath);
        assert_eq!(lonh, expected_lonh);
    }

    #[test]
    fn test_shape_to_geometry_closed_ring() {
        let shape_coord = &[
            [1.0, 2.0], 
            [3.0, 4.0], 
            [5.0, 6.0], 
            [1.0, 2.0] 
        ];

        let proj = Proj::new_known_crs("EPSG:4326", "EPSG:3857", None).unwrap();
        let polygon = shape_to_geometry(shape_coord, &proj).unwrap();

        assert!(polygon.exterior().is_closed());
        assert_eq!(polygon.exterior().coords_count(), 4);
    }

    #[test]
    fn test_shape_to_geometry_unclosed_ring() {
        let shape_coord = &[
            [1.0, 2.0], 
            [3.0, 4.0], 
            [5.0, 6.0], 
        ];

        let proj = Proj::new_known_crs("EPSG:4326", "EPSG:3857", None).unwrap();
        let polygon = shape_to_geometry(shape_coord, &proj).unwrap();
        assert_eq!(polygon.exterior().coords_count(), 4);
    }

    #[test]
    fn test_process_lat_lon() {
        // Test process lat lon with default arguments
        let path_to_nc = Path::new("example/input_ERA5/era5-crop.nc");
        let path_to_nc_str = path_to_nc.to_str().unwrap();
        let path_to_hru = Path::new("example/maps/HRUs_coarse.shp");
        let path_to_hru_str = path_to_hru.to_str().unwrap();
        let args = vec![
            "grid_gravitas",
            "-n", path_to_nc_str,
            "-d", "longitude,latitude",
            "-v", "longitude,latitude",
            "-s", "custom.shp",
            "-c", "custom_id",
            "-o", path_to_hru_str,
        ];
        let cli = Cli::parse_from(args);
        let (lath, lonh, nlat, nlon) = process_lat_lon(&cli).unwrap();
        assert_eq!(nlat, 81);
        assert_eq!(nlon, 121);
        assert_eq!(lath.shape(), &[82, 122]);
        assert_eq!(lonh.shape(), &[82, 122]);

        // Test process lat lon with grd_bnds set to true
        let args = vec![
            "grid_gravitas",
            "-n", path_to_nc_str,
            "-d", "longitude,latitude",
            "-v", "longitude,latitude",
            "-s", "custom.shp",
            "-c", "custom_id",
            "-o", path_to_hru_str,
            "--grd-bnds",
        ];
        let cli = Cli::parse_from(args);
        let (lath, lonh, nlat, nlon) = process_lat_lon(&cli).unwrap();
        assert_eq!(nlat, 80);
        assert_eq!(nlon, 120);
        assert_eq!(lath.shape(), &[81, 121]);
        assert_eq!(lonh.shape(), &[81, 121]);
    }

    #[test]
    fn test_create_grid_cells() {
        let lath: Array2<f32> = Array2::from_shape_vec((2, 2), vec![10.0, 20.0, 30.0, 40.0]).unwrap();
        let lonh: Array2<f32> = Array2::from_shape_vec((2, 2), vec![100.0, 110.0, 120.0, 130.0]).unwrap();
        let target_epsg = "EPSG:3857".to_string();
        let grid_cell_geom = create_grid_cells(1, 1, &lath, &lonh, &target_epsg).unwrap();
        assert_eq!(grid_cell_geom.len(), 1);
        assert_eq!(grid_cell_geom[0].len(), 1);
    }

    #[test]
    fn test_intersection_with_polygon_overlap() {
        // Test when there is overlap between the grid cell and the shape
        // Create a Polygon without interior rings
        let grid_cell: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let shape_polygon: Polygon<f64> = polygon![
            (x: 5.0, y: 5.0),
            (x: 5.0, y: 15.0),
            (x: 15.0, y: 15.0),
            (x: 15.0, y: 5.0),
            (x: 5.0, y: 5.0),
        ];

        // Cast Polygon into Geometry to satisfy the calculate_intersection function
        let shape_geometry = Geometry::Polygon(shape_polygon);
        let intersection = calculate_intersection(&grid_cell, &shape_geometry);
        assert!(intersection.is_some());

        // Destructure the intersection to get the MultiPolygon
        // Probably better to check equality of Polygons, but couldn't find it
        let expected_polygon = polygon![
            (x: 5.0, y: 5.0),
            (x: 5.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 5.0),
            (x: 5.0, y: 5.0),
        ];
        match intersection {
            Some(Geometry::Polygon(ref polygon)) => {
                assert!(polygon.contains(&expected_polygon));
            },
            Some(Geometry::MultiPolygon(ref mpoly)) => {
                assert_eq!(mpoly.0.len(), 1);
                let intersected_polygon = &mpoly.0[0];
                assert!(intersected_polygon.contains(&expected_polygon));
            },
            _ => panic!("Expected a Polygon or MultiPolygon geometry"),
        }
    }

    #[test]
    fn test_intersection_with_polygon_nooverlap() {
        // Test when there is no overlap between the grid cell and the shape
        let grid_cell: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let shape_polygon: Polygon<f64> = polygon![
            (x: 15.0, y: 15.0),
            (x: 15.0, y: 20.0),
            (x: 20.0, y: 20.0),
            (x: 20.0, y: 15.0),
            (x: 15.0, y: 15.0),
        ];
        let shape_geometry = Geometry::Polygon(shape_polygon);
        let intersection = calculate_intersection(&grid_cell, &shape_geometry);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_intersection_with_polygon_touching_edge() {
        // Test when we touch corners of the grid cell
        let grid_cell: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let shape_polygon: Polygon<f64> = polygon![
            (x: 10.0, y: 0.0),
            (x: 10.0, y: 10.0),
            (x: 20.0, y: 10.0),
            (x: 20.0, y: 0.0),
            (x: 10.0, y: 0.0),
        ];
        let shape_geometry = Geometry::Polygon(shape_polygon);
        let intersection = calculate_intersection(&grid_cell, &shape_geometry);
        assert!(intersection.is_some());
        let expected_linestring = LineString::from(vec![
            Coord { x: 10.0, y: 10.0 },
            Coord { x: 10.0, y: 0.0 },
        ]);
        match intersection {
            Some(Geometry::LineString(ref linestring)) => {
                assert_eq!(linestring, &expected_linestring);
            }
            _ => {
                panic!("Expected a LineString geometry");
            }
        }
    }

    #[test]
    fn test_intersection_with_non_polygon_geometry() {

        let grid_cell: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let line = geo::LineString::from(vec![
            geo::Coord { x: 5.0, y: 5.0 },
            geo::Coord { x: 15.0, y: 15.0 },
        ]);
        let line_geometry = Geometry::LineString(line);
        let intersection = calculate_intersection(&grid_cell, &line_geometry);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_intersection_with_multipolygon_overlap() {
        let grid_cell: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let shape_polygon1: Polygon<f64> = polygon![
            (x: 5.0, y: 5.0),
            (x: 5.0, y: 15.0),
            (x: 15.0, y: 15.0),
            (x: 15.0, y: 5.0),
            (x: 5.0, y: 5.0),
        ];
        let shape_polygon2: Polygon<f64> = polygon![
            (x: 15.0, y: 5.0),
            (x: 15.0, y: 15.0),
            (x: 25.0, y: 15.0),
            (x: 25.0, y: 5.0),
            (x: 15.0, y: 5.0),
        ];
        let shape_multipolygon = Geometry::MultiPolygon(MultiPolygon(vec![shape_polygon1, shape_polygon2]));
        let intersection = calculate_intersection(&grid_cell, &shape_multipolygon);
        assert!(intersection.is_some());
        let expected_polygon = polygon![
            (x: 5.0, y: 5.0),
            (x: 5.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 5.0),
            (x: 5.0, y: 5.0),
        ];
        match intersection {
            Some(Geometry::Polygon(ref polygon)) => {
                assert!(polygon.contains(&expected_polygon));
            },
            Some(Geometry::MultiPolygon(ref mpoly)) => {
                assert_eq!(mpoly.0.len(), 1);
                let intersected_polygon = &mpoly.0[0];
                assert!(intersected_polygon.contains(&expected_polygon));
            },
            _ => panic!("Expected a Polygon or MultiPolygon geometry"),
        }
    }

    #[test]
    fn test_intersection_with_multipolygon_no_overlap() {
        let grid_cell: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let shape_polygon1: Polygon<f64> = polygon![
            (x: 15.0, y: 15.0),
            (x: 15.0, y: 20.0),
            (x: 20.0, y: 20.0),
            (x: 20.0, y: 15.0),
            (x: 15.0, y: 15.0),
        ];
        let shape_polygon2: Polygon<f64> = polygon![
            (x: 25.0, y: 15.0),
            (x: 25.0, y: 20.0),
            (x: 30.0, y: 20.0),
            (x: 30.0, y: 15.0),
            (x: 25.0, y: 15.0),
        ];
        let shape_multipolygon = Geometry::MultiPolygon(MultiPolygon(vec![shape_polygon1, shape_polygon2]));
        let intersection = calculate_intersection(&grid_cell, &shape_multipolygon);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_calculate_area_with_polygon() {
        let polygon: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let geometry = Geometry::Polygon(polygon);
        let area = calculate_area(&geometry);
        assert_eq!(area, 100.0);
    }

    #[test]
    fn test_calculate_area_with_multipolygon() {
        let polygon1: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 10.0),
            (x: 10.0, y: 10.0),
            (x: 10.0, y: 0.0),
            (x: 0.0, y: 0.0),
        ];
        let polygon2: Polygon<f64> = polygon![
            (x: 10.0, y: 0.0),
            (x: 10.0, y: 10.0),
            (x: 20.0, y: 10.0),
            (x: 20.0, y: 0.0),
            (x: 10.0, y: 0.0),
        ];
        let multipolygon = MultiPolygon(vec![polygon1, polygon2]);
        let geometry = Geometry::MultiPolygon(multipolygon);
        let area = calculate_area(&geometry);
        assert_eq!(area, 200.0);
    }

    #[test]
    fn test_calculate_area_with_non_polygon() {
        let line = geo::LineString::from(vec![
            geo::Coord { x: 5.0, y: 5.0 },
            geo::Coord { x: 15.0, y: 15.0 },
        ]);
        let line_geometry = Geometry::LineString(line);
        let area = calculate_area(&line_geometry);
        assert_eq!(area, 0.0);
    }

    #[test]
    fn test_single_polygon_intersection_fully_contained() {
        // Basin polygon simple square
        let shape_polygon: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];
        let shape_geometry = Geometry::Polygon(shape_polygon.clone());
        // Single grid cell that fully contains the basin
        let grid_cell_polygon: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 4.0),
            (x: 4.0, y: 4.0),
            (x: 4.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];
        let grid_cell = GridCell {
            ilat: 0,
            ilon: 0,
            polygon: grid_cell_polygon.clone(),
        };
        
        // Building required RTree
        let rtree = RTree::bulk_load(vec![grid_cell]);
        let basin_id = String::from("test_basin");
        let nlon = 1;
        let area_error_threshold = 0.01;
        let index = 0;

        let netcdf_data_chunk = process_single_shape_netcdf(
            &shape_geometry,
            &rtree,
            nlon,
            area_error_threshold,
            index,
        );
        let txt_data_chunk = process_single_shape_txt(
            &basin_id,
            &shape_geometry,
            &rtree,
            nlon,
            area_error_threshold,
        );
        assert!(!netcdf_data_chunk.is_empty(), "NetCDF data chunk should not be empty");
        assert!(!txt_data_chunk.is_empty(), "TXT data chunk should not be empty");

        // Check expected fraction. This should be fully contained in the grid cell
        let expected_fraction = 1.0;
        assert_eq!(netcdf_data_chunk[0].0, expected_fraction);
        assert_eq!(txt_data_chunk[0].2, expected_fraction);

    }

    #[test]
    fn test_single_polygon_intersection_partial_contained() {
        // Test with two grid cells each containing half of the basin
        let shape_polygon: Polygon<f64> = polygon![
            (x: 1.0, y: 0.0),
            (x: 1.0, y: 2.0),
            (x: 3.0, y: 2.0),
            (x: 3.0, y: 0.0),
            (x: 1.0, y: 0.0)
        ];
        let shape_geometry = Geometry::Polygon(shape_polygon.clone());
        // Two grid cells 
        let grid_cell_polygon1: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];
        let grid_cell_polygon2: Polygon<f64> = polygon![
            (x: 2.0, y: 0.0),
            (x: 2.0, y: 2.0),
            (x: 4.0, y: 2.0),
            (x: 4.0, y: 0.0),
            (x: 2.0, y: 0.0)
        ];
        let grid_cell1 = GridCell {
            ilat: 0,
            ilon: 0,
            polygon: grid_cell_polygon1.clone(),
        };
        let grid_cell2 = GridCell {
            ilat: 0,
            ilon: 1,
            polygon: grid_cell_polygon2.clone(),
        };
        let rtree = RTree::bulk_load(vec![grid_cell1, grid_cell2]);
        let basin_id = String::from("test_basin");
        let nlon = 2;
        let area_error_threshold = 0.01;
        let index = 0;
        let netcdf_data_chunk = process_single_shape_netcdf(
            &shape_geometry,
            &rtree,
            nlon,
            area_error_threshold,
            index,
        );
        let txt_data_chunk = process_single_shape_txt(
            &basin_id,
            &shape_geometry,
            &rtree,
            nlon,
            area_error_threshold,
        );
        assert!(!netcdf_data_chunk.is_empty(), "NetCDF data chunk should not be empty");
        assert!(!txt_data_chunk.is_empty(), "TXT data chunk should not be empty");

        // Check expected fraction for each grid cell
        let expected_fraction = 0.5;
        assert_eq!(netcdf_data_chunk[0].0, expected_fraction);
        assert_eq!(netcdf_data_chunk[1].0, expected_fraction);
        assert_eq!(txt_data_chunk[0].2, expected_fraction);
        assert_eq!(txt_data_chunk[1].2, expected_fraction);
    }

    #[test]
    fn test_process_shape_intersections() {
        // Test with two grid cells each containing half of the basin
        let shape_polygon: Polygon<f64> = polygon![
            (x: 1.0, y: 0.0),
            (x: 1.0, y: 2.0),
            (x: 3.0, y: 2.0),
            (x: 3.0, y: 0.0),
            (x: 1.0, y: 0.0)
        ];
        let shape_geometry = Geometry::Polygon(shape_polygon.clone());
        // Two grid cells 
        let grid_cell_polygon1: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];
        let grid_cell_polygon2: Polygon<f64> = polygon![
            (x: 2.0, y: 0.0),
            (x: 2.0, y: 2.0),
            (x: 4.0, y: 2.0),
            (x: 4.0, y: 0.0),
            (x: 2.0, y: 0.0)
        ];
        let grid_cell_geom = vec![vec![grid_cell_polygon1.clone(), grid_cell_polygon2.clone()]];
        let shapes = vec![("test_basin".to_string(), shape_geometry)];
        let nlat = 1;
        let nlon = 2;
        let netcdf_data = process_shape_intersections_netcdf(nlon, grid_cell_geom.clone(), shapes.clone()).unwrap();
        let rvn_data = process_shape_intersections_txt(nlat, nlon, grid_cell_geom.clone(), shapes.clone()).unwrap();
        let expected_fraction = 0.5;
        assert_eq!(netcdf_data[0].0, expected_fraction);
        assert_eq!(netcdf_data[1].0, expected_fraction);
        assert_eq!(rvn_data.txt_data[0].2, expected_fraction);
        assert_eq!(rvn_data.txt_data[1].2, expected_fraction);
    }

    #[test]
    fn test_parallel_process_shape_intersections() {
        // Test with two grid cells each containing half of the basin
        let shape_polygon: Polygon<f64> = polygon![
            (x: 1.0, y: 0.0),
            (x: 1.0, y: 2.0),
            (x: 3.0, y: 2.0),
            (x: 3.0, y: 0.0),
            (x: 1.0, y: 0.0)
        ];
        let shape_geometry = Geometry::Polygon(shape_polygon.clone());
        // Two grid cells 
        let grid_cell_polygon1: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];
        let grid_cell_polygon2: Polygon<f64> = polygon![
            (x: 2.0, y: 0.0),
            (x: 2.0, y: 2.0),
            (x: 4.0, y: 2.0),
            (x: 4.0, y: 0.0),
            (x: 2.0, y: 0.0)
        ];
        let grid_cell_geom = vec![vec![grid_cell_polygon1.clone(), grid_cell_polygon2.clone()]];
        let shapes = vec![("test_basin".to_string(), shape_geometry)];
        let nlat = 1;
        let nlon = 2;
        let netcdf_data = parallel_process_shape_intersections_netcdf(nlon, grid_cell_geom.clone(), shapes.clone()).unwrap();
        let rvn_data = parallel_process_shape_intersections_txt(nlat, nlon, grid_cell_geom.clone(), shapes.clone()).unwrap();
        let expected_fraction = 0.5;
        assert_eq!(netcdf_data[0].0, expected_fraction);
        assert_eq!(netcdf_data[1].0, expected_fraction);
        assert_eq!(rvn_data.txt_data[0].2, expected_fraction);
        assert_eq!(rvn_data.txt_data[1].2, expected_fraction);
    }

}