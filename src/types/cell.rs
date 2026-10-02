use nalgebra::{Matrix4, Vector3};

/// Unit cell parameters: lengths in Å, angles in degrees
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Cell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl Cell {
    /// Matrix that converts fractional coordinates to cartesian coordinates
    pub fn conversion_matrix(&self) -> Matrix4<f64> {
        let Cell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } = *self;
        let cos_alpha = alpha.to_radians().cos();
        let cos_beta = beta.to_radians().cos();
        let cos_gamma = gamma.to_radians().cos();
        let sin_gamma = gamma.to_radians().sin();

        Matrix4::new(
            a,
            b * cos_gamma,
            c * cos_beta,
            0.0,
            0.0,
            b * sin_gamma,
            c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
            0.0,
            0.0,
            0.0,
            c * ((1.0 - cos_alpha.powi(2) - cos_beta.powi(2) - cos_gamma.powi(2)
                + 2.0 * cos_alpha * cos_beta * cos_gamma)
                .sqrt())
                / sin_gamma,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        )
    }

    /// Converts a fractional position to a cartesian position
    pub fn fractional_to_cartesian(&self, position: [f64; 3]) -> [f64; 3] {
        self.conversion_matrix()
            .transform_vector(&Vector3::from(position))
            .into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[test]
    fn test_orthogonal_cell() {
        let cell = Cell {
            a: 10.0,
            b: 20.0,
            c: 30.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        let position = cell.fractional_to_cartesian([0.1, 0.2, 0.3]);

        assert!(relative_eq!(position[..], [1.0, 4.0, 9.0][..], epsilon = 1.0e-9));
    }
}
