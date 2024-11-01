/*
A no_std library for gaussian curve coefficents calculation.
Copyright (C) 2024  joker2770

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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
