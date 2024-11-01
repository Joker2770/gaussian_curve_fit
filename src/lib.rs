#![no_std]

pub mod gaussian_curve_fit;

#[cfg(test)]
mod tests {
    use gaussian_curve_fit::gaussian_curve::GaussianCoefficents;

    use super::*;

    #[test]
    fn it_works() {
        let mut gaussian_coes = GaussianCoefficents::default();
        let xdata = [
            -8.0f32, -6.0f32, -4.0f32, -2.0f32, 0.0f32, 2.0f32, 4.0f32, 6.0f32,
        ];
        let ydata = [
            6.7f32, 10.6f32, 13.5f32, 15.7f32, 16.6f32, 15.4f32, 14.2f32, 10.3f32,
        ];
        let (x_arr, y_arr) = GaussianCoefficents::get_matrix_data_from_8_points(&xdata, &ydata);
        gaussian_coes.get_coefficents_from_8_matrix_data(&x_arr, &y_arr);

        assert!((gaussian_coes.value(-8.0f32) - 6.7f32).abs() < 1.0f32);
        assert!((gaussian_coes.value(-4.0f32) - 13.5f32).abs() < 1.0f32);
        assert!((gaussian_coes.value(0.0f32) - 16.6f32).abs() < 1.0f32);
        assert!((gaussian_coes.value(4.0f32) - 14.2f32).abs() < 1.0f32);
    }
}
