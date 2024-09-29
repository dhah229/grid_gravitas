use ndarray::Array2;

/// Area threshold to determine the error of the area
pub const AREA_THRESHOLD: f64 = 0.05;

pub enum OutputDataType {
    NetCDF(Vec<(f64, usize, usize)>),
    Txt(RvnGridWeights),
}

/// Struct to hold the Raven grid weights
#[derive(Debug)]
pub struct RvnGridWeights {
    pub txt_data: Vec<(String, usize, f64)>,
    pub nsubbasins: i32,
    pub ncells: i32,
}

/// Function to compare dimension names
pub fn dims_match(dims: &[String], expected: &[String]) -> bool {
    dims == expected
}

pub fn dims_match_reversed(dims: &[String], expected: &[String]) -> bool {
    dims == expected.iter().cloned().rev().collect::<Vec<String>>().as_slice()
}
/// Check if latitudes are in [-90, 90]
pub fn is_valid_lat(lat: &Array2<f32>) -> Result<(), Box<dyn std::error::Error>> {
    let is_valid_lat = lat.iter().all(|&x| x >= -90.0 && x <= 90.0);
    if !is_valid_lat {
        return Err("Latitudes must be in [-90, 90]".into());
    }
    Ok(())
}
/// Check if longitudes are in [-180, 180] or [0, 360]
pub fn is_valid_lon(lon: &Array2<f32>) -> Result<(), Box<dyn std::error::Error>> {
    let is_180_to_180 = lon.iter().all(|&x| x >= -180.0 && x <= 180.0);
    let is_0_to_360 = lon.iter().all(|&x| x >= 0.0 && x <= 360.0);
    if !(is_180_to_180 || is_0_to_360) {
        return Err("Longitudes must be in either [-180, 180] or [0, 360]".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_valid_lat() {
        let lat = Array2::from_shape_vec((2, 2), vec![0.0, 90.0, -90.0, 0.0]).unwrap();
        assert!(is_valid_lat(&lat).is_ok());
    }

    #[test]
    fn test_is_valid_lat_invalid() {
        let lat = Array2::from_shape_vec((2, 2), vec![0.0, 90.0, -90.1, 0.0]).unwrap();
        assert!(is_valid_lat(&lat).is_err());
    }

    #[test]
    fn test_is_valid_lon() {
        let lon = Array2::from_shape_vec((2, 2), vec![0.0, 180.0, -180.0, 0.0]).unwrap();
        assert!(is_valid_lon(&lon).is_ok());
    }

    #[test]
    fn test_is_valid_lon_invalid() {
        let lon = Array2::from_shape_vec((2, 2), vec![0.0, 180.0, -180.1, 0.0]).unwrap();
        assert!(is_valid_lon(&lon).is_err());
    }

    #[test]
    fn test_dims_match_identical() {
        let dims = vec!["x".to_string(), "y".to_string()];
        let expected = vec!["x".to_string(), "y".to_string()];
        assert!(dims_match(&dims, &expected));
    }

    #[test]
    fn test_dims_match_different() {
        let dims = vec!["x".to_string(), "y".to_string()];
        let expected = vec!["y".to_string(), "x".to_string()];
        assert!(!dims_match(&dims, &expected));
    }

    #[test]
    fn test_dims_match_reversed_valid() {
        let dims = vec!["x".to_string(), "y".to_string()];
        let expected = vec!["y".to_string(), "x".to_string()];
        assert!(dims_match_reversed(&dims, &expected));
    }

    #[test]
    fn test_dims_match_reversed_invalid() {
        let dims = vec!["x".to_string(), "y".to_string()];
        let expected = vec!["x".to_string(), "z".to_string()];
        assert!(!dims_match_reversed(&dims, &expected));
    }
}