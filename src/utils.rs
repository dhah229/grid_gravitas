
pub struct RvnGridWeights {
    pub txt_data: Vec<(i32, usize, f64)>,
    pub nsubbasins: i32,
    pub ncells: i32,
}

// Function to compare dimension names
pub fn dims_match(dims: &[String], expected: &[String]) -> bool {
    dims == expected
}

pub fn dims_match_reversed(dims: &[String], expected: &[String]) -> bool {
    dims == expected.iter().cloned().rev().collect::<Vec<String>>().as_slice()
}

#[cfg(test)]
mod tests {
    use super::*;

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