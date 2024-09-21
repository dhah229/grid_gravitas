mod io;
mod geometry;
mod utils;
mod cli;

use std::{
    error::Error,
    path::Path,
};

use clap::Parser;

use io::{read_shapefile, write_netcdf_output, write_txt_output};
use geometry::{
    process_lat_lon, 
    create_grid_cells, 
    process_shape_intersections,
};
use cli::Cli;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    // Make grid of latitudes and longitudes
    let (lath, lonh, nlat, nlon) = process_lat_lon(&args)?;

    // Create grid of polygons for each grid cell
    let grid_cell_geom = create_grid_cells(nlat, nlon, &lath, &lonh)?;

    // Get the shapefile ID and geometry to a HashMap
    let shapes = read_shapefile(Path::new(&args.shp), &args.col, "3573")?;

    // Process the shape intersections with grid cells
    let (netcdf_data, mut rvn_data) = process_shape_intersections(
        nlat, nlon, 
        &grid_cell_geom, 
        &shapes,
    )?;

    // Write to output file
    let output_path = Path::new(&args.out);
    if args.rv_out {
        write_txt_output(output_path, &mut rvn_data)?;
    } else {
        write_netcdf_output(output_path, &netcdf_data)?;
    }
    Ok(())
}


