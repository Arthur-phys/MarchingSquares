mod lib;
use crate::lib::MarchingSquares;


fn main() {
    fn identityb(b: (f32,f32)) -> f32 {
        b.0
    }
    let a = MarchingSquares::new(40.,5.,0.5,identityb);
    println!("{:?}",a.values_matrix);
    println!("{:?}",a.boolean_matrix);
    println!("{:?}",a.index_matrix);
    println!("{:?}",a.interpol_matrix);
    a.jsondump("../Plotting/final.txt","Campo eléctrico de un pez");
}