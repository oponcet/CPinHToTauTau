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
        light_blue="#bbd3f3",
        bright_blue="#0aa4f6",
        dark_blue="#10376c",
        light_purple="#cbabfa",
        bright_purple="#bf5ffa",
        purple="#761fc3",
        aubergine="#431f67",
        yellow="#ebe341",
        dark_yellow="#b1b111",
        bright_orange="#e96a0c",
        dark_orange="#8c4106",
        red="#FF0000",
        teal="#008080",
        grey="#808080",
        maroon="#800000",
        brown="#8f3007",
        black="#000000",
    )

    if (p := config.get_process("h_ggf_tautau", default=None)):
        p.color1 = cfg.x.colors.black
        p.label = (
            r"$H_{ggf} \rightarrow \tau\tau$"
        )
        
    if (p := config.get_process("h", default=None)):
        p.color1 = cfg.x.colors.purple

    if (p := config.get_process("tt", default=None)):
        p.color1 = cfg.x.colors.light_purple
        p.label = r"$t\bar{t}$"

    if (p := config.get_process("st", default=None)):
        p.color1 = cfg.x.colors.aubergine

    #if (p := config.get_process("dy", default=None)):
    #    p.color1 = cfg.x.colors.brown
    if (p := config.get_process("dy_lep_m50", default=None)):
        p.color1 = cfg.x.colors.brown
        p.label = r"$Z\to \ell \ell (e/\mu/\tau)$"
    if (p := config.get_process("dy_z2ll", default=None)):
        p.color1 = cfg.x.colors.light_blue
        p.label = r"$Z\to \ell \ell ~(ee/\mu\mu)$"
    if (p := config.get_process("dy_z2tautau", default=None)):
        p.color1 = cfg.x.colors.dark_blue
        p.label = r"$Z\to \tau \tau$"

    if (p := config.get_process("dy_lep_m10to50", default=None)):
        p.color1 = cfg.x.colors.grey
        p.label = r"$Z\to \ell \ell (M < 50)$"

    if (p := config.get_process("vv", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("vvv", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("multiboson", default=None)):
        p.color1 = cfg.x.colors.yellow

    if (p := config.get_process("w", default=None)):
        p.color1 = cfg.x.colors.bright_orange
        p.label = "W"
    if (p := config.get_process("w_lnu", default=None)):
        p.color1 = cfg.x.colors.bright_orange
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
