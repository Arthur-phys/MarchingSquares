use ndarray::{Array2,Array,Ix2};
use std::fs::File;
use std::io::prelude::*;

pub struct MarchingSquares{
    square_lenght: f32, // El cuadrado se genera desde el centro y hacia los lados con esta longitúd
    partition_size: f32, // El tamaño de la partición que se aplicará sobre el cuadrado con 2square_lenght de longitúd
    isovalue: f32, // El valor de la curva isoterma que queremos crear
    pub values_matrix: Array<f32,Ix2>, // La matríz con los valores de la función de la que queremos crear la isoterma
    pub boolean_matrix: Array<u8,Ix2>, // La matríz con los valores binarios que ayudan a crear la isoterma
    pub index_matrix: Array<u8,Ix2>, // La matríz con los valores del 0-15 para crear los casos de las isotermas
    pub interpol_matrix: Array<Vec<(f32,f32)>,Ix2>,
    pub mathematical_function: fn((f32,f32)) -> f32,
}

impl MarchingSquares {

    pub fn new(square_lenght: f32, partition_size: f32, isovalue: f32, mathematical_function: fn((f32,f32)) -> f32 ) -> MarchingSquares {

        let maf = |tup: (usize,usize)| { // A la función 'mathematical_function' le ponemos el formato que pide 'from_shape_fn' que es de tomar una tupla de usize y regresar un flotante 
            let floatpoint = MarchingSquares::center_point(partition_size,square_lenght,tup.0,tup.1);
            mathematical_function(floatpoint)
        };
    
        let values_matrix = Array2::from_shape_fn((partition_size as usize,partition_size as usize),maf); // Esta función resuleve un montón de mis problemas
        let boolean_matrix = MarchingSquares::boolean_matrix(&values_matrix, partition_size, isovalue);
        let index_matrix = MarchingSquares::index_matrix(partition_size,&boolean_matrix);
        let interpol_matrix = MarchingSquares::interpolation_matrix(square_lenght,partition_size,&index_matrix,isovalue,mathematical_function);
        
        MarchingSquares {
            square_lenght,
            partition_size,
            isovalue,
            values_matrix,
            boolean_matrix,
            index_matrix,
            interpol_matrix,
            mathematical_function,
        }
    }

    fn index_matrix(partition_size: f32,boolean_matrix: &Array<u8,Ix2>) -> Array<u8,Ix2> { // La matríz de índices del 0 al 15

        let mut im = Array2::<u8>::zeros((partition_size as usize - 1,partition_size as usize - 1)); // im = index matrix
        
        for i in 0..(partition_size as usize - 1){
            for j in 0..(partition_size as usize - 1){
                let i1 = boolean_matrix[[i,j]]*(2_u8.pow(3));
                let i2 = boolean_matrix[[i,j+1]]*(2_u8.pow(2));
                let i3 = boolean_matrix[[i+1,j+1]]*2_u8;
                let i4 = boolean_matrix[[i+1,j]];
                im[[i,j]] = i1 + i2 + i3 + i4;
            }
        }
        im
    }

    pub fn interpolation_matrix(lenght: f32, partition_size: f32, index_matrix: &Array<u8,Ix2>, isovalue:f32, maf: fn((f32,f32)) -> f32) -> Array<Vec<(f32,f32)>,Ix2> { // Se encarga de los casos y la interpolacion
        
        let mut im: Array<Vec<(f32,f32)>,Ix2> = Array2::default((partition_size as usize - 1,partition_size as usize - 1)); // im = interpolation matrix
        
        for i in 0..(partition_size as usize - 1) {
            for j in 0..(partition_size as usize - 1) {
                if index_matrix[(i,j)] == 0 { // Nada Adentro, no regresa mas que una tupla de ceros
                    im[(i,j)] = vec![(0.,0.)]; // Los elementos de matríz vecinos están constriñidos por los nodos que comparten. En este caso no regresamos nada porque los vecinos se encargan de todo el trabajo
                } else if index_matrix[(i,j)] == 1 { // Punto inferior izquierdo dentro

                    let inside_point = MarchingSquares::center_point(partition_size,lenght,i+1,j); // El punto que queda adentro
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j); // punto auxiliar superior izquierdo (tiene la misma coordenada y que inside_point) // Estos puntos ayudan a hacer la interpolacion
                    let a = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // Punto auxiliar inferior derecho (tiene la misma coordenada x que inside_point)

                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point.1,inside_point.0,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(a.0,a.1,inside_point.0,inside_point.1,isovalue,maf);

                    im[(i,j)] = vec![(inside_point.0,qy),(px,inside_point.1)] // Los puntos tienen orden de izquierda a derecha

                } else if index_matrix[(i,j)] == 2 { // Punto inferior derecho dentro

                    let inside_point = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // El punto que queda adentro
                    let b = MarchingSquares::center_point(partition_size,lenght,i+1,j); // punto auxiliar inferior izquierdo
                    let a = MarchingSquares::center_point(partition_size,lenght,i,j+1); // Punto auxiliar superior derecho

                    let qx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point.0,inside_point.1,isovalue,maf);
                    let py = MarchingSquares::linear_interpolation_y(a.1,a.0,inside_point.1,inside_point.0,isovalue,maf);

                    im[(i,j)] = vec![(qx,inside_point.1),(inside_point.0,py)];

                } else if index_matrix[(i,j)] == 3 { // Puntos inferiores dentro

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Los puntos dentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i+1,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j); // Los puntos superiores que no están dentro
                    let a = MarchingSquares::center_point(partition_size,lenght,i,j+1); 
                    
                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let py = MarchingSquares::linear_interpolation_y(a.1,a.0,inside_point2.1,inside_point2.0,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(inside_point2.0,py)];

                } else if index_matrix[(i,j)] == 4 { // Punto superior derecho dentro

                    let inside_point = MarchingSquares::center_point(partition_size,lenght,i,j+1); // EL punto dentro
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j); // Punto auxiliar superior izquierdo
                    let a = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // Punto auziliar inferior derecho

                    let qx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point.0,inside_point.1,isovalue,maf); 
                    let py = MarchingSquares::linear_interpolation_y(a.1,a.0,inside_point.1,inside_point.0,isovalue,maf); 

                    im[(i,j)] = vec![(qx,inside_point.1),(inside_point.0,py)];

                } else if index_matrix[(i,j)] == 5 { // Puntos inferior izquierdo y superior derecho dentro

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Los puntos dentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j); // Punto auxiliar superior izquierdo
                    let a = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // Punto auxiliar inferior derecho

                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point2.0,inside_point2.1,isovalue,maf);
                    let rx = MarchingSquares::linear_interpolation_x(a.0,a.1,inside_point1.0,inside_point1.1,isovalue,maf);
                    let sy = MarchingSquares::linear_interpolation_y(a.1,a.0,inside_point2.1,inside_point2.0,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(px,inside_point2.1),(rx,inside_point1.1),(inside_point2.0,sy)];

                } else if index_matrix[(i,j)] == 6 { // Puntos derechos dentro

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j+1); // Los puntos dentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i+1,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j); // Punto auxiliar superior izquierdo
                    let a = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Punto auxiliar inferior izquierdo

                    let qx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point1.0,inside_point1.1,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(a.0,a.1,inside_point2.0,inside_point2.1,isovalue,maf);

                    im[(i,j)] = vec![(qx,inside_point1.1),(px,inside_point2.1)]; // Los puntos se ordenan de arriba para abajo si no se puede de izquierda a derecha

                } else if index_matrix[(i,j)] == 7 { // Todos los puntos dentro excepto por el superior izquierdo

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Los puntos dentro: inferior izquierdo y superior derecho (falta uno pero ese no ayuda para esto)
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j); // Punto auxiliar superior izquierdo
                    
                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point2.0,inside_point2.1,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(px,inside_point2.1)];

                } else if index_matrix[(i,j)] == 8 { // Punto superior izquierdo dentro (estos casos son todos simétricos a los anteriores y por eso haré copy paste de casi todos sin cambiar comentarios)

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j);
                    let b = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Puntos auxiliares inferior izquierdo y superior derecho
                    let a = MarchingSquares::center_point(partition_size,lenght,i,j+1);
                    
                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(a.0,a.1,inside_point1.0,inside_point1.1,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(px,inside_point1.1)];

                } else if index_matrix[(i,j)] == 9 { // Puntos izquierdos dentro

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j); // Los puntos dentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i+1,j);
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j+1); // Punto auxiliar superior derecho
                    let a = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // Punto auxiliar inferior derecho

                    let qx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point1.0,inside_point1.1,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(a.0,a.1,inside_point2.0,inside_point2.1,isovalue,maf);

                    im[(i,j)] = vec![(qx,inside_point1.1),(px,inside_point2.1)]; // Los puntos se ponen de arriba para abajo si no se puede de izquierda a derecha

                } else if index_matrix[(i,j)] == 10 { //Puntos superior izquierdo e inferior derecho dentro

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j); // Los puntos dentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i+1,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Punto auxiliar superior izquierdo
                    let a = MarchingSquares::center_point(partition_size,lenght,i,j+1); // Punto auxiliar inferior derecho

                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(a.0,a.1,inside_point1.0,inside_point1.1,isovalue,maf);
                    let rx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point2.0,inside_point2.1,isovalue,maf);
                    let sy = MarchingSquares::linear_interpolation_y(a.1,a.0,inside_point2.1,inside_point2.0,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(px,inside_point1.1),(rx,inside_point2.1),(inside_point2.0,sy)];


                } else if index_matrix[(i,j)] == 11 { // Todos dentro excepto por el superior derecho

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j); // Puntos dentro que ayudan (no son todos justo como en el caso 4)
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i+1,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i,j+1); // Punto auxiliar superior derecho

                    let qx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point1.0,inside_point1.1,isovalue,maf); 
                    let py = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point2.1,inside_point2.0,isovalue,maf); 

                    im[(i,j)] = vec![(qx,inside_point1.1),(inside_point2.0,py)];

                } else if index_matrix[(i,j)] == 12 { // Puntos superiores dentro

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j); // Los puntos dentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i,j+1);
                    let b = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Los puntos inferiores que no están dentro
                    let a = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); 
                    
                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let py = MarchingSquares::linear_interpolation_y(a.1,a.0,inside_point2.1,inside_point2.0,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(inside_point2.0,py)];

                } else if index_matrix[(i,j)] == 13 { // Todos los puntos dentro excepto por el inferior derecho 

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i+1,j); // El punto que queda adentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i,j+1); // punto auxiliar inferior izquierdo
                    let b = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // Punto auxiliar superior derecho

                    let qx = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point1.0,inside_point1.1,isovalue,maf);
                    let py = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point2.1,inside_point2.0,isovalue,maf);

                    im[(i,j)] = vec![(qx,inside_point1.1),(inside_point2.0,py)];

                } else if index_matrix[(i,j)] == 14 { // Todos los puntos dentro excepto por el inferior izquierdo

                    let inside_point1 = MarchingSquares::center_point(partition_size,lenght,i,j); // El punto que queda adentro
                    let inside_point2 = MarchingSquares::center_point(partition_size,lenght,i+1,j+1); // punto auxiliar superior izquierdo (tiene la misma coordenada y que inside_point) // Estos puntos ayudan a hacer la interpolacion
                    let b = MarchingSquares::center_point(partition_size,lenght,i+1,j); // Punto auxiliar inferior derecho (tiene la misma coordenada x que inside_point)

                    let qy = MarchingSquares::linear_interpolation_y(b.1,b.0,inside_point1.1,inside_point1.0,isovalue,maf);
                    let px = MarchingSquares::linear_interpolation_x(b.0,b.1,inside_point2.0,inside_point2.1,isovalue,maf);

                    im[(i,j)] = vec![(inside_point1.0,qy),(px,inside_point2.1)] // Los puntos tienen orden de izquierda a derecha

                } else if index_matrix[(i,j)] == 15 { // Da las esquinas de la celda de la matriz (Todo queda adentro)
                    im[(i,j)] = vec![(0.,0.)];
                }
            }
        }
        im
    }

    fn boolean_matrix(values_matrix: &Array<f32,Ix2>,partition_size: f32, isovalue: f32) -> Array<u8,Ix2> {
        
        let mut bm = Array2::<u8>::zeros((partition_size as usize, partition_size as usize));

        for index in values_matrix.indexed_iter() {
            if values_matrix[index.0] > isovalue {
                bm[index.0] = 1;
            } else {
                bm[index.0] = 0;
            }
        }

        bm
    }

    pub fn center_point(psize: f32, lenght: f32, i: usize, j: usize) -> (f32,f32) { // Regresa el centro de una celda de la matríz para usar sobre nuestra función de dos variables
        (lenght*(1.-((2*j) as f32 + 1.)/psize),lenght*(((2*i) as f32 + 1.)/psize - 1.)) // psize = partition size // i === x, j === y. PPodría haber un error aquí porque supuse mal el orden de i y j
    }

    fn linear_interpolation_x(a:f32,b:f32,c:f32,d:f32,isovalue:f32,mathfuncborrowed: fn((f32,f32)) -> f32) -> f32 {// Acá va la definición de la interpolación lineal para usar en todos lados, pero requiere de una función, por lo que a lo mejor es buena idea poner la función como variable miembro en MarchingSquares
        a + (c-a)*((isovalue-mathfuncborrowed((a,b)))/(mathfuncborrowed((c,d))-mathfuncborrowed((a,b))))
    } 
    fn linear_interpolation_y(a:f32,b:f32,c:f32,d:f32,isovalue:f32,mathfuncborrowed: fn((f32,f32)) -> f32) -> f32 { // Debe haber dos funciones de interpolacion porque los argumentos de f(x,y) siempre van en ese orden
        a + (c-a)*((isovalue-mathfuncborrowed((b,a)))/(mathfuncborrowed((d,c))-mathfuncborrowed((b,a))))
    }

    pub fn jsondump(&self,filename:&str,title:&str) {
        let mut vector_from_matrix: Vec<Vec<(f32,f32)>> = vec![];
        for index in self.interpol_matrix.indexed_iter() {
            vector_from_matrix.push(index.1.to_vec());
        }
        let string_from_vector = format!("{{\"points\":{:?},\n\"length\":{},\n\"title\":\"{}\"}}",vector_from_matrix,self.square_lenght,title);
        let mut file = File::create(filename).unwrap();
        file.write(string_from_vector.as_bytes()).unwrap();
    }
}