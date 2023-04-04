#include <iostream>
#include "CLHEP/Random/Random.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/Random.h"

int main() {
	  // Inicializar el generador de números aleatorios
	     CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	  //
	  //     // Generar un número aleatorio uniforme entre 0 y 1
	         double x = CLHEP::HepRandom::getFlat();
	  //
	  //         // Calcular la raíz cuadrada del número aleatorio
	             double y = std::sqrt(x);
	  //
	  //             // Imprimir el número aleatorio y su raíz cuadrada
	                 std::cout << "Número aleatorio: " << x << std::endl;
	                   std::cout << "Raíz cuadrada: " << y << std::endl;
	  //
	                     return 0;
	                     }
	  //
