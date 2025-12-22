# Grid Gravitas
`grid_gravitas` is a high-performance Rust tool designed to map geospatial points or shapes to grid cells in a NetCDF file. It leverages R-Trees for fast spatial lookups and supports parallel processing to handle large datasets efficiently.


## Background
In hydrological modelling, it is often necessary to compute the spatial averages of gridded data using catchment (or basin) areas. To achieve this, one of the most common ways is to pre-compute the weights of the grid cells that intersect with the basin polygons and apply those weights to determine the spatial averages. This is particularly useful when you have a large dataset (e.g., across multiple decades and multiple variables) and you want to avoid unnecessary computation. 

For small number of polygons and/or small grids, I recommend the [GridWeightsGenerator](https://github.com/julemai/GridWeightsGenerator) or [xESMF](https://xesmf.readthedocs.io/en/stable/), which are simple to use and easy to get started. However, for numerous polygons and/or large grids, grid weights generation can take a long time and can be quite burdensome if we have to compute the weights for multiple grids (e.g., ERA5, HRDPS, HRRR, and so on). 


## Installation

To install this package, it is assumed that you have [Rust](https://www.rust-lang.org/tools/install) installed on your machine, along with `GDAL`, `GEOS`, `PROJ` and `netcdf` binaries. The package can be installed through docker (please see below). 

To install the package, you can run the following cargo command:
```bash
cargo build --release
```
This will create the binary in the `target/release` folder. You can copy this binary to a folder in your `PATH` or run it directly from the `target/release` folder.

### Docker (optional)
You can build and run a containerized version using the GDAL base image:
```bash
# Build the image
docker build -t grid-gravitas .

# Run with local data mounted at /data (example inputs provided in example/)
docker run --rm -v $(pwd):/data grid-gravitas \
  --nc /data/example/input_ERA5/era5-crop.nc \
  --shp /data/example/maps/HRUs_coarse.shp \
  --out /data/output.nc
```
The container entrypoint is `grid_gravitas`, and it will print the help text by default.

## Usage
See [example.txt](example.txt) for an example. Note the example dataset comes from the [GridWeightsGenerator](https://github.com/julemai/GridWeightsGenerator).

Here is the help message for the tool:
```
Processes NetCDF and shapefiles

Usage: grid_gravitas [OPTIONS] --nc <NC> --shp <SHP>

Options:
  -n, --nc <NC>        Path to the NetCDF file
  -d, --dim <DIM>      Dimension names of longitude (x) and latitude (y) (in this order). Example: "rlon,rlat", or "x,y" [default: x y]
  -c, --coord <COORD>  Variable names of longitude and latitude in NetCDF (in this order). Example: "lon,lat" [default: lon lat]
  -s, --shp <SHP>      Path to the shapefile
  -i, --id <ID>        Name of the id in the shapefile to use as the key [default: ID]
  -b, --grd-bnds       Flag if coordinates refer to grid centers (default) or grid bounds
  -o, --out <OUT>      Path to the output file [default: output.nc]
  -r, --rv-out         Flag for Raven grid weights file
  -p, --parallel       Flag for parallel processing
  -e, --epsg <EPSG>    Target EPSG to reproject the grid and shapes for area calculation [default: EPSG:8857]
  -v, --verbose        Verbose mode
  -h, --help           Print help
  -V, --version        Print version
```
In general, the tool requires a NetCDF file that has the coordinates of the grid, and a shapefile that has the polygons that you want to generate the weights for. For example, if you have a NetCDF file of [HRDPS](https://eccc-msc.github.io/open-data/msc-data/nwp_hrdps/readme_hrdps-datamart_en/) data (size (2540,1290) with rotated-pole (rlon, rlat) as dimensions), it requires that the data also has the longitude and latitude coordinate variables (lon, lat) to calculate the weights (see [example.txt](example.txt) for data in rotated-pole grid). By default, a NetCDF file is created that stores the grid weights in the same format as xESMF (see [to_netcdf](https://xesmf.readthedocs.io/en/stable/user_api.html#xesmf.frontend.Regridder.to_netcdf)). Additionally, to use the weights in Raven, you can use the `-r` flag to output the weights in the Raven format. It should be noted that there are example datasets in the [example](example/) folder, which was retrieved from the [GridWeightsGenerator](https://github.com/julemai/GridWeightsGenerator).

For HRDPS data with ~7000 basins across Canada, the tool takes about 10 minutes to calculate the weights in parallel on a 16-core machine. 

```
grid_gravitas --nc HRDPS_grid.nc --shp basins.shp --dimname rlon,rlat --varname lon,lat --id station_no --out HRDPS_weights.nc --epsg EPSG:3573 --parallel
```

## Docker Support

Using Docker is the recommended way to develop and run `grid_gravitas` as it bundles all necessary geospatial C-libraries (GDAL, PROJ, GEOS, and NetCDF).

### Running Tests
To run the internal Rust test suite and verify the environment:
```bash
docker-compose up test

### Running the Tool
```
docker run --rm -v "$(pwd):/data" grid_gravitas \
    --nc /data/grid.nc \
    --shp /data/basins.shp \
    --id station_id \
    --out /data/weights.nc \
    --parallel
```
