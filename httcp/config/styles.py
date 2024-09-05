# coding: utf-8

"""
Style definitions.
"""

import order as od

from columnflow.util import DotDict


def stylize_processes(config: od.Config) -> None:
    """
    Adds process colors and adjust labels.
    """
    cfg = config

    # recommended cms colors
    cfg.x.colors = DotDict(
        bright_blue="#3f90da",
        dark_blue="#011c87",
        purple="#832db6",
        aubergine="#964a8b",
        yellow="#f7c331",
        bright_orange="#ffa90e",
        dark_orange="#e76300",
        red="#bd1f01",
        teal="#92dadd",
        grey="#94a4a2",
        brown="#a96b59",
    )

    if (p := config.get_process("h_ggf_tautau", default=None)):
        p.color1 = cfg.x.colors.dark_blue
        p.label = (
            r"$H_{ggf} \rightarrow \tau\tau$"
        )
        
    if (p := config.get_process("h", default=None)):
        p.color1 = cfg.x.colors.purple

    if (p := config.get_process("tt", default=None)):
        p.color1 = cfg.x.colors.bright_orange
        p.label = r"$t\bar{t}$"

    if (p := config.get_process("st", default=None)):
        p.color1 = cfg.x.colors.aubergine

    if (p := config.get_process("dy", default=None)):
        p.color1 = cfg.x.colors.dark_orange
    if (p := config.get_process("dy_lep_m50", default=None)):
        p.color1 = cfg.x.colors.dark_orange
        p.label = r"$Z\to \ell \ell$"
    if (p := config.get_process("dy_lep_m10to50", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$Z\to \ell \ell (LM)$"

    if (p := config.get_process("vv", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("vvv", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("multiboson", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("w", default=None)):
        p.color1 = cfg.x.colors.teal
        p.label = "W"
    if (p := config.get_process("w_lnu", default=None)):
        p.color1 = cfg.x.colors.teal
        p.label = "W"

    if (p := config.get_process("ewk", default=None)):
        p.color1 = cfg.x.colors.brown

    if (p := config.get_process("ttv", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$t\bar{t} + V$"

    if (p := config.get_process("ttvv", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$t\bar{t} + VV$"

    if (p := config.get_process("tt_multiboson", default=None)):
        p.color1 = cfg.x.colors.grey

    if (p := config.get_process("qcd", default=None)):
        p.color1 = cfg.x.colors.red
