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
    pub dimname: Vec<String>,

    /// Variable names of longitude and latitude in NetCDF (in this order). Example: "lon,lat"
    #[arg(short = 'v', long, default_values = &["lon", "lat"], value_delimiter = ',')]
    pub varname: Vec<String>,

    /// Path to the shapefile
    #[arg(short = 's', long)]
    pub shp: String,

    /// Name of the column in the shapefile to use as the key
    #[arg(short = 'c', long, default_value = "ID")]
    pub col: String,

    /// Path to the output file
    #[arg(short = 'o', long, default_value = "output.nc")]
    pub out: String,

    /// Flag for Raven grid weights file
    #[arg(short = 'r', long)]
    pub rv_out: bool,

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_values() {
        let args = vec!["grid_gravitas", "-n", "file.nc", "-s", "file.shp"];
        let cli = Cli::parse_from(args);

        assert_eq!(cli.nc, "file.nc");
        assert_eq!(cli.dimname, vec!["x", "y"]);
        assert_eq!(cli.varname, vec!["lon", "lat"]);
        assert_eq!(cli.shp, "file.shp");
        assert_eq!(cli.col, "ID");
        assert_eq!(cli.out, "output.nc");
        assert_eq!(cli.rv_out, false);
    }

    #[test]
    fn test_custom_values() {
        let args = vec![
            "grid_gravitas",
            "-n", "custom.nc",
            "-d", "rlon,rlat",
            "-v", "longitude,latitude",
            "-s", "custom.shp",
            "-c", "custom_id",
            "-o", "custom_output.nc",
            "-r",
        ];
        let cli = Cli::parse_from(args);

        assert_eq!(cli.nc, "custom.nc");
        assert_eq!(cli.dimname, vec!["rlon", "rlat"]);
        assert_eq!(cli.varname, vec!["longitude", "latitude"]);
        assert_eq!(cli.shp, "custom.shp");
        assert_eq!(cli.col, "custom_id");
        assert_eq!(cli.out, "custom_output.nc");
        assert_eq!(cli.rv_out, true);
    }
}