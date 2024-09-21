use ndarray::{s, Array, Array2, Zip};
use proj::Proj;

use std::{
    error::Error,
    path::Path,
    collections::HashMap,
};

use geo::{
    algorithm::area::Area,
    BooleanOps,
    BoundingRect,
    Geometry,
    HasDimensions,
    Intersects,
    MultiPolygon,
    Coord,
    LineString,
    MapCoords,
    Polygon,
};

use crate::utils::RvnGridWeights;
use crate::io::read_lat_lon;
use crate::cli::Cli;



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

    let mut exterior = LineString::from(coords);

    // Transform the coordinates if necessary
    exterior = exterior
        .try_map_coords(|coord| proj.convert(coord)).expect("Failed to transform coordinates");

    let polygon = Polygon::new(exterior, vec![]);

    Ok(polygon)
}


/// Process lat/lon and create grid cells
pub fn process_lat_lon(args: &Cli) -> Result<(Array2<f32>, Array2<f32>, usize, usize), Box<dyn Error>> {
    let file = netcdf::open(Path::new(&args.nc)).expect("Failed to open file.");
    let (lat_2d, lon_2d) = read_lat_lon(&file, &args).expect("Failed to read lat/lon");

    let (lath, lonh) = create_gridcells_from_centers(&lat_2d, &lon_2d);
    
    let nlat = lath.nrows() - 1;
    let nlon = lonh.ncols() - 1;

    Ok((lath, lonh, nlat, nlon))
}

/// Create grid cells and store them as geometries
pub fn create_grid_cells(
    nlat: usize, 
    nlon: usize, 
    lath: &Array2<f32>, 
    lonh: &Array2<f32>
) -> Result<Vec<Vec<Polygon<f64>>>, Box<dyn Error>> {
    let mut grid_cell_geom: Vec<Vec<Polygon<f64>>> = Vec::with_capacity(nlat);
    let proj = Proj::new_known_crs("EPSG:4326", "EPSG:3573", None).expect("Failed to create projection");

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

/// Calculate the intersection between the grid cell and the shape geometry
fn calculate_intersection(grid_cell: &Polygon<f64>, shape_geometry: &Geometry) -> Option<Geometry> {
    match shape_geometry {
        Geometry::Polygon(ref shape_poly) => {
            let intersection = grid_cell.intersection(shape_poly);
            if !intersection.is_empty() {
                Some(Geometry::MultiPolygon(intersection))
            } else {
                None
            }
        },
        Geometry::MultiPolygon(ref shape_mpoly) => {
            let grid_cell_mpoly = MultiPolygon(vec![grid_cell.clone()]);
            let intersection = grid_cell_mpoly.intersection(shape_mpoly);
            if !intersection.is_empty() {
                Some(Geometry::MultiPolygon(intersection))
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
        _ => 0.0,
    }
}


/// Process the shape intersections with grid cells and collect results
pub fn process_shape_intersections(
    nlat: usize, 
    nlon: usize, 
    grid_cell_geom: Vec<Vec<Polygon<f64>>>, 
    shapes: HashMap<String, geo_types::Geometry>, 
) -> Result<(Vec<(f64, usize, usize)>, RvnGridWeights), Box<dyn Error>> {
    let nsubbasins = shapes.len() as i32;
    let ncells = (nlat * nlon) as i32;
    let mut netcdf_data = Vec::new();
    let mut txt_data = Vec::new();

    // 5 percent area error threshold
    let area_error_threshold = 0.05;

    for (index, (basin_id, shape_geometry)) in shapes.into_iter().enumerate() {
        let shape_area = match shape_geometry {
            Geometry::Polygon(ref poly) => poly.unsigned_area(),
            Geometry::MultiPolygon(ref mpoly) => mpoly.unsigned_area(),
            _ => continue,
        };

        let shape_envelope = shape_geometry.bounding_rect().ok_or("Failed to get bounding rectangle")?;
        let mut area_all = 0.0;
        let mut cellid_fraction = Vec::new();

        for ilat in 0..nlat {
            for ilon in 0..nlon {
                let grid_cell = &grid_cell_geom[ilat][ilon];
                let grid_envelope = grid_cell.bounding_rect().ok_or("Failed to get bounding rectangle")?;

                if grid_envelope.intersects(&shape_envelope) {
                    let intersection = calculate_intersection(grid_cell, &shape_geometry);
                    if let Some(intersection_geom) = intersection {
                        let area_intersect = calculate_area(&intersection_geom);
                        if area_intersect > 0.0 {
                            area_all += area_intersect;
                            let cell_id = ilat * nlon + ilon;
                            let fraction = area_intersect / shape_area;
                            cellid_fraction.push((cell_id, fraction));
                        }
                    }
                }
            }
        }
        // Calculate error
        let error = (shape_area - area_all) / shape_area;

        // Cloning basin_id shouldn't be too bad since they should be short.
        if error.abs() > area_error_threshold && shape_area > 500_000.0 {
            println!("Error for basin {}: {:.2}%", basin_id, error * 100.0);
            for (cell_id, fraction) in cellid_fraction {
                netcdf_data.push((fraction, (cell_id + 1), (index + 1)));
                txt_data.push((basin_id.clone(), cell_id, fraction))
            }
        } else {
            // Adjust such that weights sum up to 1.0
            let correction_factor = 1.0 / (1.0 - error);
            for (cell_id, fraction) in cellid_fraction {
                let corrected = fraction * correction_factor;
                netcdf_data.push((corrected, (cell_id + 1), (index + 1)));
                txt_data.push((basin_id.clone(), cell_id, corrected))
            }
        }
    }
    let rvn_data = RvnGridWeights {
        txt_data,
        nsubbasins,
        ncells,
    };

    Ok((netcdf_data, rvn_data))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use geo::CoordsIter; 

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

}