mod io;
mod geometry;
mod utils;
mod cli;
mod rtree;

use std::{
    error::Error,
    path::Path,
};

use clap::Parser;

use io::{read_shapefile, write_netcdf_output, write_txt_output};
use geometry::{
    process_lat_lon, 
    create_grid_cells, 
    process_shape_intersections_netcdf,
    process_shape_intersections_txt,
    parallel_process_shape_intersections_netcdf,
    parallel_process_shape_intersections_txt,
};
use cli::Cli;
use utils::OutputDataType;


fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    run_with_args(args)
}

// Why 8857?
// We need a projected coordinate system to calculate areas and 8857 is an Equal Earth Projection.
// Although this has inaccuracies near the poles, it is good enough for our purposes.
// We could use something like ESRI:102017(North Pole Lambert Azimuthal Equal Area). 
// But doesn't cover the globe.

fn run_with_args(args: Cli) -> Result<(), Box<dyn Error>> {
    // Make grid of latitudes and longitudes
    let (lath, lonh, nlat, nlon) = process_lat_lon(&args)?;

    // Create grid of polygons for each grid cell
    let grid_cell_geom = create_grid_cells(nlat, nlon, &lath, &lonh)?;

    // Get the shapefile ID and geometry to a Vec of tuples
    let shapes = read_shapefile(Path::new(&args.shp), &args.col, "8857")?;

    // Process the shape intersections with grid cells
    let data: OutputDataType = match (args.parallel, args.rv_out) {
        (true, true) => OutputDataType::Txt(parallel_process_shape_intersections_txt(
            nlat, nlon, grid_cell_geom, shapes,
        )?),
        (true, false) => OutputDataType::NetCDF(parallel_process_shape_intersections_netcdf(
            nlon, grid_cell_geom, shapes,
        )?),
        (false, true) => OutputDataType::Txt(process_shape_intersections_txt(
            nlat, nlon, grid_cell_geom, shapes,
        )?),
        (false, false) => OutputDataType::NetCDF(process_shape_intersections_netcdf(
            nlon, grid_cell_geom, shapes,
        )?),
    };

    // Write to output file
    let output_path = Path::new(&args.out);
    match data {
        OutputDataType::NetCDF(data) => {
            write_netcdf_output(output_path, &data)?;
        }
        OutputDataType::Txt(mut data) => {
            write_txt_output(output_path, &mut data)?;
        }
    }
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use std::fs::File;
    use std::io::BufRead;

    #[test]
    fn test_e2e_series_nc() {
        // Test E2E
        let temp_file = NamedTempFile::new().expect("Failed to create temporary file");
        let temp_path = temp_file.path().to_path_buf();
        let temp_path_str = temp_path.to_str().unwrap();
        let path_to_nc = Path::new("example/input_ERA5/era5-crop.nc");
        let path_to_nc_str = path_to_nc.to_str().unwrap();
        let path_to_hru = Path::new("example/maps/HRUs_coarse.shp");
        let path_to_hru_str = path_to_hru.to_str().unwrap();
        let args = Cli {
            nc: path_to_nc_str.to_string(),
            dimname: vec!["longitude".to_string(), "latitude".to_string()],
            varname: vec!["longitude".to_string(), "latitude".to_string()],
            shp: path_to_hru_str.to_string(),
            col: "HRU_ID".to_string(),
            grd_bnds: false,
            out: temp_path_str.to_string(),
            rv_out: false,
            parallel: false,
        };
        run_with_args(args).expect("Failed to run with args");

        let nc_file = netcdf::open(temp_path).expect("Failed to open netcdf file");
        let dim = nc_file.dimension("n_s").expect("Dimension not found");
        assert_eq!(dim.len(), 123);
    }

    #[test]
    fn test_e2e_parallel_nc() {
        // Test E2E
        let temp_file = NamedTempFile::new().expect("Failed to create temporary file");
        let temp_path = temp_file.path().to_path_buf();
        let temp_path_str = temp_path.to_str().unwrap();
        let path_to_nc = Path::new("example/input_ERA5/era5-crop.nc");
        let path_to_nc_str = path_to_nc.to_str().unwrap();
        let path_to_hru = Path::new("example/maps/HRUs_coarse.shp");
        let path_to_hru_str = path_to_hru.to_str().unwrap();
        let args = Cli {
            nc: path_to_nc_str.to_string(),
            dimname: vec!["longitude".to_string(), "latitude".to_string()],
            varname: vec!["longitude".to_string(), "latitude".to_string()],
            shp: path_to_hru_str.to_string(),
            col: "HRU_ID".to_string(),
            grd_bnds: false,
            out: temp_path_str.to_string(),
            rv_out: false,
            parallel: true,
        };
        run_with_args(args).expect("Failed to run with args");

        let nc_file = netcdf::open(temp_path).expect("Failed to open netcdf file");
        let dim = nc_file.dimension("n_s").expect("Dimension not found");
        assert_eq!(dim.len(), 123);
    }


    #[test]
    fn test_e2e_series_rv() {
        // Test E2E
        let temp_file = NamedTempFile::new().expect("Failed to create temporary file");
        let temp_path = temp_file.path().to_path_buf();
        let temp_path_str = temp_path.to_str().unwrap();
        let path_to_nc = Path::new("example/input_ERA5/era5-crop.nc");
        let path_to_nc_str = path_to_nc.to_str().unwrap();
        let path_to_hru = Path::new("example/maps/HRUs_coarse.shp");
        let path_to_hru_str = path_to_hru.to_str().unwrap();
        let args = Cli {
            nc: path_to_nc_str.to_string(),
            dimname: vec!["longitude".to_string(), "latitude".to_string()],
            varname: vec!["longitude".to_string(), "latitude".to_string()],
            shp: path_to_hru_str.to_string(),
            col: "HRU_ID".to_string(),
            grd_bnds: false,
            out: temp_path_str.to_string(),
            rv_out: true,
            parallel: false,
        };
        run_with_args(args).expect("Failed to run with args");

        let output_file = File::open(&temp_path).expect("Failed to open temporary file");
        let reader = std::io::BufReader::new(output_file);
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
        let lines_trimmed: Vec<String> = lines.iter().map(|l| l.trim().to_string()).collect();

        // Check that the number of HRUs is 51 and the number of grid cells is 9801
        let hru_line = lines_trimmed.iter().find(|line| line.contains("NumberHRUs") && line.contains("51"));
        let grid_cells_line = lines_trimmed.iter().find(|line| line.contains("NumberGridCells") && line.contains("9801"));

        assert!(hru_line.is_some(), "NumberHRUs line not found or value does not match");
        assert!(grid_cells_line.is_some(), "NumberGridCells line not found or value does not match");


    }

    #[test]
    fn test_e2e_parallel_rv() {
        // Test E2E
        let temp_file = NamedTempFile::new().expect("Failed to create temporary file");
        let temp_path = temp_file.path().to_path_buf();
        let temp_path_str = temp_path.to_str().unwrap();
        let path_to_nc = Path::new("example/input_ERA5/era5-crop.nc");
        let path_to_nc_str = path_to_nc.to_str().unwrap();
        let path_to_hru = Path::new("example/maps/HRUs_coarse.shp");
        let path_to_hru_str = path_to_hru.to_str().unwrap();
        let args = Cli {
            nc: path_to_nc_str.to_string(),
            dimname: vec!["longitude".to_string(), "latitude".to_string()],
            varname: vec!["longitude".to_string(), "latitude".to_string()],
            shp: path_to_hru_str.to_string(),
            col: "HRU_ID".to_string(),
            grd_bnds: false,
            out: temp_path_str.to_string(),
            rv_out: true,
            parallel: true,
        };
        run_with_args(args).expect("Failed to run with args");

        let output_file = File::open(&temp_path).expect("Failed to open temporary file");
        let reader = std::io::BufReader::new(output_file);
        let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
        let lines_trimmed: Vec<String> = lines.iter().map(|l| l.trim().to_string()).collect();

        // Check that the number of HRUs is 51 and the number of grid cells is 9801
        let hru_line = lines_trimmed.iter().find(|line| line.contains("NumberHRUs") && line.contains("51"));
        let grid_cells_line = lines_trimmed.iter().find(|line| line.contains("NumberGridCells") && line.contains("9801"));

        assert!(hru_line.is_some(), "NumberHRUs line not found or value does not match");
        assert!(grid_cells_line.is_some(), "NumberGridCells line not found or value does not match");


    }
    
}

