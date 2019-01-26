#include "DataProcessUsingPython.hh"
#include <iostream>
#include <Python.h>
#include <vector>

//This function is largely copied from python tumor webstie: 
//https://docs.python.org/2/extending/embedding.html

void DataProcessUsingPython::PlotCellSurvivalCurve(string fileName, string color)
{
    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    int i;
    
    std::vector<std::string>  vec_str;
    vec_str.push_back("plotCellSF");
    vec_str.push_back("plotCellSF");
    vec_str.push_back("plotSFCurve");

    vec_str.push_back(fileName.c_str()); // the input argument for PyString_FromString is const char*, so use this line
    vec_str.push_back(color.c_str());
//     vec_str.push_back("'cellState.csv'");
//     vec_str.push_back("'b'");
//     vec_str.push_back("multiply");
//     vec_str.push_back("multiply");
//     vec_str.push_back("multiply");
//     vec_str.push_back("9");
//     vec_str.push_back("3");
    char** argv1  = (char**)&vec_str[0];
    int argc = vec_str.size();
    cerr<<"RADCellSimulationInitializePyWrapper argv  [0]="<<argv1[0]<<endl;
    
    for (int i = 0 ; i < vec_str.size() ; ++i)
    {
        cerr<<"argument["<<i<<"]="<<argv1[i]<<endl;
    }

    if (argc < 3) 
    {
        fprintf(stderr,"Usage: call pythonfile funcname [args]\n");

    }

    Py_Initialize();
    pName = PyString_FromString(argv1[1]);
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv1[2]);
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc - 3);
            for (i = 0; i < argc - 3; ++i) {
//                 pValue = PyString_FromString(atoi(argv1[i + 3]));// if need number argument, use this line
                pValue = PyString_FromString((argv1[i + 3])); // if need string argument, use this line
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                printf("Result of call: %ld\n", PyInt_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv1[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv1[1]);
    }
    Py_Finalize();

}

