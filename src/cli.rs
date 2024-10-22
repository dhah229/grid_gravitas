use clap::Parser;

#[derive(Parser)]
#[command(name = "grid_gravitas")]
#[command(author = "David Hah")]
#[command(version = "1.0")]
#[command(about = "Processes NetCDF and shapefiles", long_about = None)]
pub struct Cli {
    /// Path to the NetCDF file
    #[arg(short = 'n', long)]
    pub nc: String,

    /// Dimension names of longitude (x) and latitude (y) (in this order). Example: "rlon,rlat", or "x,y"
    #[arg(short = 'd', long, default_values = &["x", "y"], value_delimiter = ',')]
    pub dim: Vec<String>,

    /// Variable names of longitude and latitude in NetCDF (in this order). Example: "lon,lat"
    #[arg(short = 'c', long, default_values = &["lon", "lat"], value_delimiter = ',')]
    pub coord: Vec<String>,

    /// Path to the shapefile
    #[arg(short = 's', long)]
    pub shp: String,

    /// Name of the id in the shapefile to use as the key
    #[arg(short = 'i', long, default_value = "ID")]
    pub id: String,

    /// Flag if coordinates refer to grid centers (default) or grid bounds
    #[arg(short = 'b', long)]
    pub grd_bnds: bool,

    /// Path to the output file
    #[arg(short = 'o', long, default_value = "output.nc")]
    pub out: String,

    /// Flag for Raven grid weights file
    #[arg(short = 'r', long)]
    pub rv_out: bool,

    /// Flag for parallel processing
    #[arg(short = 'p', long)]
    pub parallel: bool,

    /// Target EPSG to reproject the grid and shapes for area calculation
    #[arg(short = 'e', long, default_value = "EPSG:8857")]
    pub epsg: String,

    /// Verbose mode
    #[arg(short = 'v', long)]
    pub verbose: bool,

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_values() {
        let args = vec!["grid_gravitas", "-n", "file.nc", "-s", "file.shp"];
        let cli = Cli::parse_from(args);

        assert_eq!(cli.nc, "file.nc");
        assert_eq!(cli.dim, vec!["x", "y"]);
        assert_eq!(cli.coord, vec!["lon", "lat"]);
        assert_eq!(cli.shp, "file.shp");
        assert_eq!(cli.id, "ID");
        assert_eq!(cli.out, "output.nc");
        assert_eq!(cli.rv_out, false);
    }

    #[test]
    fn test_custom_values() {
        let args = vec![
            "grid_gravitas",
            "-n", "custom.nc",
            "-d", "rlon,rlat",
            "-c", "longitude,latitude",
            "-s", "custom.shp",
            "-i", "custom_id",
            "-o", "custom_output.nc",
            "-r",
            "-p",
            "-v",
        ];
        let cli = Cli::parse_from(args);

        assert_eq!(cli.nc, "custom.nc");
        assert_eq!(cli.dim, vec!["rlon", "rlat"]);
        assert_eq!(cli.coord, vec!["longitude", "latitude"]);
        assert_eq!(cli.shp, "custom.shp");
        assert_eq!(cli.id, "custom_id");
        assert_eq!(cli.out, "custom_output.nc");
        assert_eq!(cli.rv_out, true);
        assert_eq!(cli.parallel, true);
        assert_eq!(cli.verbose, true);
    }
}