cdef extern from "MeshTypes.h":
    cdef cppclass ElementClass:
        ElementClass()  except +
        ElementClass(int dim_) except +

    int getElement(ElementClass &element)
