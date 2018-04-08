from typing import *
import logging
import numpy as np

try:
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri

    def convert_r_obj(v: Any, obj_to_obj: bool=True, verbose: bool=True) -> Any:
        """Function with manually specified conversion from a r-object to a python object
        """
        if type(v) == ro.rinterface.RNULLType:
            return None
        elif type(v) == ro.vectors.Matrix:
            return np.array(v)
        elif type(v) == ro.vectors.FloatVector:
            return np.array(v, dtype="float64")
        elif type(v) == ro.vectors.IntVector:
            return np.array(v, dtype="int64")
        elif type(v) == ro.rinterface.RNULLType:
            return None
        elif type(v) == ro.vectors.ListVector:
            try:
                return {v.names[i]: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
            except TypeError:
                return {i: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
        elif type(v) == ro.vectors.StrVector:
            if len(v) == 1:
                return str(v[0])
            else:
                try:
                    return {v.names[i]: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
                except TypeError:
                    return {i: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
        elif type(v) == ro.vectors.DataFrame:
            from rpy2.robjects import pandas2ri
            return pandas2ri.ri2py(v)
        elif type(v) == ro.methods.RS4:
            if obj_to_obj:
                class RS4Object(object):
                    def __repr__(self) -> str:
                        return f"< RS4Object with attributes: {list(self.__dict__.keys())} >"
                rs4obj = RS4Object()
                for k in tuple(v.slotnames()):
                    setattr(rs4obj, k, convert_r_obj(v.slots[k], obj_to_obj=obj_to_obj))
                return rs4obj
            else:
                return {k: convert_r_obj(v.slots[k]) for k in tuple(v.slotnames())}
        else:
            if type(v) != str:
                if verbose:
                    print(f"not supported yet {type(v)}")
            return v

except:
    pass
