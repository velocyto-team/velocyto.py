from typing import Any
from contextlib import suppress
import numpy as np
from loguru import logger

with suppress(ImportError):
    import rpy2.robjects as ro

    def convert_r_obj(v: Any, obj_to_obj: bool = True, verbose: bool = True) -> Any:
        """Function with manually specified conversion from a r-object to a python object"""
        match type(v):
            case ro.rinterface.RNULLType:
                return None
            case ro.vectors.Matrix:
                return np.array(v)
            case ro.vectors.FloatVector:
                return np.array(v, dtype="float64")
            case ro.vectors.IntVector:
                return np.array(v, dtype="int64")
            case ro.vectors.ListVector:
                try:
                    return {v.names[i]: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
                except TypeError:
                    return {i: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
            case ro.vectors.StrVector:
                if len(v) == 1:
                    return str(v[0])
                try:
                    return {v.names[i]: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
                except TypeError:
                    return {i: convert_r_obj(v[i], obj_to_obj=obj_to_obj) for i in range(len(v))}
            case ro.vectors.DataFrame:
                from rpy2.robjects import pandas2ri

                return pandas2ri.ri2py(v)
            case ro.methods.RS4:
                if not obj_to_obj:
                    return {k: convert_r_obj(v.slots[k]) for k in tuple(v.slotnames())}

                class RS4Object(object):
                    def __repr__(self) -> str:
                        return f"< RS4Object with attributes: {list(self.__dict__.keys())} >"

                rs4obj = RS4Object()
                for k in tuple(v.slotnames()):
                    setattr(rs4obj, k, convert_r_obj(v.slots[k], obj_to_obj=obj_to_obj))
                return rs4obj
            case _:
                if type(v) != str and verbose:
                    logger.info(f"not supported yet {type(v)}")
                return v
