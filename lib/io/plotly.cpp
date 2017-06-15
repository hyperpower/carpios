#include "plotly.h"

namespace carpio{

Plotly_actor::Plotly_actor(){
	Py_Initialize();
}

Plotly_actor::~Plotly_actor(){
	//Py_Finalize();
};

}
