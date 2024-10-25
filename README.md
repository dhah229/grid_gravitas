# Grid Gravitas
Grid Gravitas is a Rust-based command line tool designed to calculate weights based on the intersection of shapefiles and a grid for spatial averages. 

## Background
In hydrological modelling, it is often necessary to compute the spatial averages of gridded data using catchment (or basin) areas. To achieve this, one of the most common ways is to pre-compute the weights of the grid cells that intersect with the basin polygons and reuse those weights to determine the spatial averages. This is particularly useful when you have a large dataset (e.g., across multiple decades and multiple variables) and you want to avoid unnecessary computation. 

For small number of polygons and/or small grids, I recommend the [GridWeightsGenerator](https://github.com/julemai/GridWeightsGenerator) or [xESMF](https://xesmf.readthedocs.io/en/stable/), which are simple to use and easy to get started. However, for numerous polygons and/or large grids, grid weights generation can take a long time and can be quite burdensome if we have to compute the weights for multiple grids (e.g., ERA5, HRDPS, HRRR, and so on). This inspired me to build this tool that takes inspiration from the aforementioned methods but leverages parallel computation targted for large grid weights generation. Under the hood, Grid Gravitas uses GDAL, GEOS, PROJ, and netcdf for data processing and computation.

## Installation
To install this package, it is assumed that you have [Rust](https://www.rust-lang.org/tools/install) installed on your machine, along with `GDAL`, `GEOS`, `PROJ` and `netcdf` binaries. For reference, you can install these dependencies with `conda` with the environment file provided in the [environment.yml](environment.yml) file. You may need to set the necessary environment variables like `LIBCLANG_PATH`, `PKG_CONFIG_PATH`, and `LD_LIBRARY_PATH` (see for example [env_vars.sh](env_vars.sh)).

To install the package, you can run the following cargo command:
```bash
cargo build --release
```
This will create the binary in the `target/release` folder. You can copy this binary to a folder in your `PATH` or run it directly from the `target/release` folder.

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


## TODO
- [ ] Dockerize the tool for easier installation
- [ ] Add more examples and documentation
- [ ] Add more tests
- [ ] Remove the nested Vec<>

