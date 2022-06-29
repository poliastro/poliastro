from enum import Flag, auto


class PropagatorKind(Flag):
    ELLIPTIC = auto()
    PARABOLIC = auto()
    HYPERBOLIC = auto()
