import json
import re
from pathlib import Path
from typing import TypedDict


class SBSVariable(TypedDict):
    """Interpreted metadata for a single SBS variable.

    :param name: converted variable name with special characters removed
    :param aliases: other names that return this variable
    :param format: ``%f`` for a numeric value, ``%s`` for a text string
    :param units: variable units; empty string if undefined
    :param title: variable title, such as Temperature or Conductivity
    :param info: other info, such as ITS-90, Secondary, or sensor model
    """

    name: str
    aliases: list[str]
    format: str
    units: str
    title: str
    info: list[str]


def _build_lookup() -> dict:
    """Build a lookup table of the variables in sbs_variables.json
    """
    definitions_path = Path(__file__).parent / "resources" / "sbs_variables.json"
    lookup: dict = {}
    definitions = json.loads(definitions_path.read_text(encoding="utf-8"))
    for entry in definitions:
        for alias in entry["aliases"]:
            lookup[alias] = entry
    return lookup


# run once on import
_LOOKUP = _build_lookup()


def interpret_sbs_variable(alias: str) -> SBSVariable:
    """SBS variable name interpretation.

    Converts the strings that define SBS variables to strings that exclude
    special characters. Gives the units, title, and sensor type for each
    variable and defines whether it is a number or a text string.

    All variables are defined in the ``sbs_variables.json`` resources file,
    originally imported from SBEDataProcessing_7.26.4.pdf. There may be
    additional variables added in future versions of Sea-Bird Data Processing
    modules that are not included there.

    :param alias: text string with the SBS variable name

    :return: a typed dict with the following fields: name, aliases,
        format, units, title, and info
    """
    sbs_var = _LOOKUP.get(alias)
    if sbs_var is None:
        # Replace non-alphanummeric characters (\W) with underscores
        converted = re.sub(r"\W", "_", alias)
        return SBSVariable(
            name=converted,
            aliases=[],
            format=r"%f",
            units="",
            title=converted,
            info=[],
        )

    return SBSVariable(
        # a null name means it's the same as the alias
        name=sbs_var["name"] if sbs_var["name"] is not None else alias,
        aliases=sbs_var["aliases"],
        format=sbs_var.get("format", r"%f"),
        units=sbs_var["units"],
        title=sbs_var["title"],
        info=sbs_var["info"],
    )
