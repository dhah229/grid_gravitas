use geo::prelude::*;
use geo_types::Polygon;
use rstar::{RTree, RTreeObject, AABB};

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