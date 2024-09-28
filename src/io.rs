use netcdf;
use geo_types::Geometry as GeoGeometry;
use ndarray::{Array2, Ix1, Ix2};
use wkt::Wkt;
use std::{
    convert::TryFrom,
    error::Error,
    path::Path,
    io::{BufWriter, Write},
    fs::File,
};

use gdal::{
    Dataset as GdalDataset,
    spatial_ref::{CoordTransform, SpatialRef},
    vector::LayerAccess,
    vector::OGRFieldType,
};

use crate::utils::{dims_match, dims_match_reversed, RvnGridWeights};
use crate::geometry::meshgrid;
use crate::cli::Cli;

// Read latitude and longitude from the NetCDF file
pub fn read_lat_lon(file: &netcdf::File, args: &Cli) -> Result<(ndarray::Array2<f32>, ndarray::Array2<f32>), Box<dyn Error>> {

    // Unpack dimension and variable names
    let (var_lon, var_lat) = (&args.varname[0], &args.varname[1]);
    
    let lat_var = file
        .variable(var_lat).ok_or("Latitude variable not found")?;
    let lon_var = file
        .variable(var_lon).ok_or("Longitude variable not found")?;

    let lat_dims = lat_var.dimensions();
    let lon_dims = lon_var.dimensions();


    // Ensure both variables have the same number of dimensions
    if lat_dims.len() != lon_dims.len() {
        return Err(format!(
            "Dimension mismatch: Latitude has {} dimensions, but Longitude has {} dimensions.",
            lat_dims.len(),
            lon_dims.len()
        ).into());
    }

    // Extract dimension names as Vec<String>
    let lat_dim_names: Vec<String> = lat_dims.iter().map(|d| d.name().to_string()).collect();
    let lon_dim_names: Vec<String> = lon_dims.iter().map(|d| d.name().to_string()).collect();

    // Extract expected dimension names from CLI
    let expected_dimname = &args.dimname;

    // Initialize variables to hold the final arrays
    let lat_final: Array2<f32>;
    let lon_final: Array2<f32>;

    // Check if the dimensions are 1D
    match lat_dims.len() {
        1 => {
            // If the dimensions are 1D, convert to 2D using meshgrid
            let lat = lat_var
                .get::<f32, _>(..)?
                .into_dimensionality::<Ix1>()
                .expect("Latitude is not 1D");
            let lon = lon_var
                .get::<f32, _>(..)?
                .into_dimensionality::<Ix1>()
                .expect("Longitude is not 1D");

            let (lat_2d, lon_2d) = meshgrid(&lat, &lon);
            lat_final = lat_2d;
            lon_final = lon_2d;
        },
        2 => {
            // If the dimensions are 2D, extract the data and transpose if necessary
            // [      1      2      3   ...     1*nlon
            //   nlon+1 nlon+2 nlon+3   ...     2*nlon
            //      ...    ...    ...   ...     ...
            //      ...    ...    ...   ...  nlat*nlon ]
            let lat_2d = lat_var
                .get::<f32, _>(..)?
                .into_dimensionality::<Ix2>()
                .expect("Latitude is not 2D");
            let lon_2d = lon_var
                .get::<f32, _>(..)?
                .into_dimensionality::<Ix2>()
                .expect("Longitude is not 2D");
            
            // Check dimension order and transpose if necessary
            if dims_match(&lon_dim_names, expected_dimname) {
                // Dimensions are in the expected order, transpose
                lon_final = lon_2d.t().to_owned();
            } else if dims_match_reversed(&lon_dim_names, expected_dimname) {
                // Dimensions are in reverse order, do nothing
                lon_final = lon_2d;
            } else {
                return Err("Dimension names do not match".into());
            }

            // Repeat the same for latitude
            if dims_match(&lat_dim_names, expected_dimname) {
                // Dimensions are in the expected order, transpose
                lat_final = lat_2d.t().to_owned();
            } else if dims_match_reversed(&lat_dim_names, expected_dimname) {
                // Dimensions are in reverse order, do nothing
                lat_final = lat_2d;
            } else {
                return Err("Dimension names do not match".into());
            }
        },
        _ => {
            return Err("Coordinate variables must have the same number of dimensions (either 1 or 2)".into());
        }
    }
    Ok((lat_final, lon_final))
}

// Reads shapefile and returns a Vec of IDs and geometries
pub fn read_shapefile(
    path: &Path,
    key_colname: &str,
    target_epsg: &str,
) -> Result<Vec<(String, GeoGeometry<f64>)>, Box<dyn Error>> {
    // Open the shapefile dataset
    let dataset = GdalDataset::open(path)?;

    // Get the layer
    let mut layer = dataset.layer(0)?;

    // Get the source spatial reference
    let source_srs = layer.spatial_ref().ok_or("Layer has no spatial reference")?;

    // Create the target spatial reference
    let target_srs_def = format!("EPSG:{}", target_epsg);
    let target_srs = SpatialRef::from_definition(&target_srs_def)?;

    // Create a coordinate transformation
    let coord_transform = CoordTransform::new(&source_srs, &target_srs)?;

    let mut geometries = Vec::new();

    // Iterate over features
    for feature_result in layer.features() {
        let feature = feature_result;

        // Get the geometry as a reference
        let geometry_ref = feature.geometry().ok_or("Feature has no geometry")?;

        // Clone the geometry to get an owned, mutable copy
        let mut geometry = geometry_ref.clone();

        // Transform the geometry to the target CRS
        geometry.transform_inplace(&coord_transform)?;

        // Convert GDAL geometry to WKT
        let wkt_str = geometry.wkt()?;

        // Parse the WKT string into a Wkt object
        let wkt_obj: Wkt<f64> = Wkt::from_str(&wkt_str)?;

        // Convert Wkt object into geo_types::Geometry
        let geo_geometry = GeoGeometry::<f64>::try_from(wkt_obj)?;

        // Get the basin ID (adjust field type and name as needed)
        let field_value = feature.field(key_colname)?.ok_or("Field not found")?;

        let basin_id = match field_value.ogr_field_type() {
            OGRFieldType::OFTString => {
                field_value.into_string().ok_or("Failed to get string value")?
            },
            _ => {
                field_value.into_int().ok_or("Failed to get integer value")?.to_string()
            },
        };
        geometries.push((basin_id, geo_geometry));
    }

    Ok(geometries)
}


/// Write the output to netCDF format
pub fn write_netcdf_output(file_path: &Path, netcdf_data: &Vec<(f64, usize, usize)>) -> Result<(), Box<dyn Error>> {
    // Unpack the tuple data into separate vectors
    let s_f64: Vec<f64> = netcdf_data.iter().map(|(s, _, _)| *s).collect();
    let col_usize: Vec<usize> = netcdf_data.iter().map(|(_, col, _)| *col).collect();
    let row_usize: Vec<usize> = netcdf_data.iter().map(|(_, _, row)| *row).collect();

    // Convert the vectors to the correct data types
    let s: Vec<f32> = s_f64.iter().map(|&x| x as f32).collect();
    let col: Vec<i32> = col_usize.iter().map(|&x| x as i32).collect();
    let row: Vec<i32> = row_usize.iter().map(|&x| x as i32).collect();

    // Create the netCDF file
    let mut nc_out = netcdf::create(file_path)?;
    nc_out.add_dimension("n_s", s.len())?;
    
    // Add the variables to the netCDF file
    let mut s_var = nc_out.add_variable::<f32>("S", &["n_s"])?;
    s_var.put_values(&s, ..)?;
    
    let mut col_var = nc_out.add_variable::<i32>("col", &["n_s"])?;
    col_var.put_values(&col, ..)?;
    
    let mut row_var = nc_out.add_variable::<i32>("row", &["n_s"])?;
    row_var.put_values(&row, ..)?;

    Ok(())
}

/// Write intersection data to a .txt file
pub fn write_txt_output(output_path: &Path, rvn_data: &mut RvnGridWeights) -> Result<(), Box<dyn Error>> {
    // Sort txt_data by the integer value of the first element (converted from String)
    rvn_data.txt_data.sort_by(|a, b| {
        // Parse the first element of each tuple to i32
        let a_id: i32 = a.0.parse().expect("Failed to parse basin_id to i32");
        let b_id: i32 = b.0.parse().expect("Failed to parse basin_id to i32");
        a_id.cmp(&b_id)
    });

    // Create the output file
    let output_file = File::create(output_path)?;
    let mut writer = BufWriter::new(output_file);
    writeln!(writer, ":GridWeights")?;
    writeln!(writer, "   #")?;
    writeln!(writer, "   # [# HRUs]")?;
    writeln!(writer, "   :NumberHRUs       {}", rvn_data.nsubbasins)?;
    writeln!(writer, "   :NumberGridCells  {}", rvn_data.ncells)?;
    writeln!(writer, "   #")?;
    writeln!(writer, "   # [HRU ID] [Cell #] [w_kl]")?;
    
    for (basin_id, cell_id, fraction) in &rvn_data.txt_data {
        writeln!(writer, "   {}   {}   {}", basin_id, cell_id, fraction)?;
    }
    writeln!(writer, ":EndGridWeights")?;

    Ok(())
}



#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_read_shapefile_valid() {
        let path = Path::new("example/maps/HRUs_coarse.shp");
        let result = read_shapefile(
           path, 
            "HRU_ID", 
            "8857",
        );
        
        assert!(result.is_ok());
        let shapefile = result.unwrap();
        assert_eq!(shapefile.len(), 51);
    }

    #[test]
    fn test_read_shapefile_nonexistent() {
        let path = Path::new("example/nonexistent.shp");
        let result = read_shapefile(
           path, 
            "HRU_ID", 
            "8857",
        );
        
        assert!(result.is_err());
    }

    #[test]
    fn test_read_lat_lon_valid() {
        let path = Path::new("example/input_ERA5/era5-crop.nc");
        let file = netcdf::open(path).unwrap();
        let args = Cli {
            nc: "example/input_ERA5/era5-crop.nc".to_string(),
            dimname: vec!["longitude".to_string(), "latitude".to_string()],
            varname: vec!["longitude".to_string(), "latitude".to_string()],
            shp: "example/maps/HRUs_coarse.shp".to_string(),
            col: "HRU_ID".to_string(),
            out: "output.nc".to_string(),
            rv_out: false,
            grd_bnds: false,
            parallel: false,
        };
        let result = read_lat_lon(&file, &args);        
        assert!(result.is_ok());
        let (lat, lon) = result.unwrap();
        assert_eq!(lat.shape(), [81,121]);
        assert_eq!(lon.shape(), [81,121]);

    }
}