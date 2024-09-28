use geo::prelude::*;
use geo_types::Polygon;
use rstar::{RTree, RTreeObject, AABB};

// These functions are used to create the grid cells and build the RTree
// The RTree helps to quickly find the grid cell that intersects with a shape instead of iterating over all grid cells

// GridCell struct to store the grid cell index and its geometry
#[derive(Debug)]
pub struct GridCell {
    pub ilat: usize,
    pub ilon: usize,
    pub polygon: Polygon<f64>,
}

impl RTreeObject for GridCell {
    type Envelope = AABB<[f64; 2]>;

    fn envelope(&self) -> Self::Envelope {
        let bbox = self.polygon.bounding_rect().expect("Grid cell must have a bounding rectangle");
        AABB::from_corners([bbox.min().x, bbox.min().y], [bbox.max().x, bbox.max().y])
    }
}

impl PartialEq for GridCell {
    fn eq(&self, other: &Self) -> bool {
        self.ilat == other.ilat && self.ilon == other.ilon
    }
}

impl Eq for GridCell {}

// Probably need to change this we are not doing vec of vecs.
pub fn build_grid_rtree(grid_cell_geom: Vec<Vec<Polygon<f64>>>) -> RTree<GridCell> {
    let mut grid_cells = Vec::new();

    for (ilat, row) in grid_cell_geom.into_iter().enumerate() {
        for (ilon, polygon) in row.into_iter().enumerate() {
            grid_cells.push(GridCell {
                ilat,
                ilon,
                polygon: polygon,
            });
        }
    }

    RTree::bulk_load(grid_cells)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::{Polygon, polygon};

    #[test]
    fn test_rtree_construction() {
        let polygon1: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 1.0),
            (x: 1.0, y: 1.0),
            (x: 1.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];
        let polygon2: Polygon<f64> = polygon![
            (x: 1.0, y: 1.0),
            (x: 1.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 1.0),
            (x: 1.0, y: 1.0)
        ];
        let grid_cell_geom: Vec<Vec<Polygon>> = vec![
            vec![polygon1],
            vec![polygon2]
        ];
        let rtree = build_grid_rtree(grid_cell_geom);
        assert_eq!(rtree.size(), 2);
    }

    #[test]
    fn test_rtree_bounding_box() {
        let polygon: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];

        let grid_cell = GridCell {
            ilat: 0,
            ilon: 0,
            polygon: polygon.clone(),
        };

        // Get the envelope (bounding box) of the grid cell
        let envelope = grid_cell.envelope();

        // The bounding box should match the polygon's extent
        assert_eq!(envelope.lower(), [0.0, 0.0]);
        assert_eq!(envelope.upper(), [2.0, 2.0]);
    }

    #[test]
    fn test_gridcell_equality() {
        let polygon1: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 1.0),
            (x: 1.0, y: 1.0),
            (x: 1.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];

        let polygon2: Polygon<f64> = polygon![
            (x: 0.0, y: 0.0),
            (x: 0.0, y: 2.0),
            (x: 2.0, y: 2.0),
            (x: 2.0, y: 0.0),
            (x: 0.0, y: 0.0)
        ];

        // Two grid cells with the same ilat and ilon but different polygons
        let grid_cell_1 = GridCell {
            ilat: 0,
            ilon: 0,
            polygon: polygon1.clone(),
        };

        let grid_cell_2 = GridCell {
            ilat: 0,
            ilon: 0,
            polygon: polygon2.clone(),
        };

        // Test equality based on ilat and ilon only (polygons are ignored in the eq implementation)
        assert_eq!(grid_cell_1, grid_cell_2);
    }

}