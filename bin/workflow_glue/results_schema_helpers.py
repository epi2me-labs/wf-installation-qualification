#!/usr/bin/env python
"""Helpers for the auto-generated schema code."""
import pandas as pd
from pydantic import BaseModel as PydanticBaseModel


class BaseModel(PydanticBaseModel):
    """Extend base model."""

    class Config:
        """Config items for the pydantic code."""

        # make enums json serializable
        use_enum_values = True

    @staticmethod
    def to_dataframe(list_of_objects):
        """Make a dataframe from a list of objects."""
        if type(list_of_objects) is not list:
            raise TypeError(f"Object must be instance of {list}")
        return pd.DataFrame.from_dict(
            [dict(obj) for obj in list_of_objects])
