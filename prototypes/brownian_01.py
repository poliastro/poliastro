#!/usr/bin/env python
# coding: utf-8

# # First Prototype, `OrbitArray` Class

# ## Objective: User-facing API
# 
# From a user's perspective, the ultimate goal is to achieve something along the following lines:
# 
# ```python
# from poliastro.twobody import Orbit, OrbitArray
# from poliastro.ephem import Ephem
# 
# from typing import Any
# 
# from astropy.time import Time, TimeDelta
# from datetime import datetime, timedelta
# 
# import numpy as np
# 
# ## Some test-data
# 
# time_zero: datetime = Time('1970-01-01 00:00').datetime
# days: int = 36525 # 100a
# times: Time = Time([Time(time_zero + timedelta(days = day)) for day in range(0, days)])
# 
# data: list[dict[str, Any]] = {...} # a list of classical orbital elements (and attractors)
# 
# ## Poliastro ingest, type 1
# 
# orbs: list[Orbit] = [
#     Orbit.from_classical(**item) # ingest data
#     for item in data
# ]
# 
# orbarr: OrbitArray = OrbitArray(orbs)
# 
# ## Poliastro ingest, type 2
# 
# orbarr: OrbitArray = OrbitArray.from_classical(
#     **{key: [item[key] for item in data] for k in data[0].keys()}
# )
# 
# ## Time-deltas, 1D array
# 
# dt1d: TimeDelta = times - orbarr.epochs
# 
# assert len(dt1d) == len(orbarr) == len(orbarr.epochs)
# 
# ## Propagate array with time-deltas
# 
# orbarr_propagated: OrbitArray = orbarr.propagate(dt1d)
# 
# assert len(orbarr_propagated) == len(orbarr)
# 
# ## Propagate array with array of time stamps (hiding the time-deltas)
# 
# # Consistent with Orbit.propagate as it returns another OrbitArray
# orbarr_propagated: OrbitArray = orbarr.propagate(times)
# 
# assert len(orbarr_propagated) == len(orbarr)
# 
# ## Time-deltas, 2D array
# 
# dt2d = np.repeat(times, len(orbarr), axis = 0) - orbarr.epochs[:, None]
# 
# assert dt2d.shape == (len(orbarr), len(times))
# 
# # Juanlu's idea: A special `propagate_array` function:
# 
# from poliastro.somewhere import propagate_array
# 
# xxx: UnkownType = propagate_array(orbarr, astro_timedeltas)
# 
# ## gh-1364: Juanlu's preferred solution, extended (longer-term idea?)
# 
# ephem: Ephem = Ephem.from_orbit(orbs[0], **kwargs)
# ephemarr: EphemArray = EphemArray.from_orbitarray(orbarr, **kwargs)
# ```

# ## Strategy
# 
# To **simplify** matters, this is not astrodynamics but **Brownian motion**, of sorts. This notebook is looking at the relationships between `Orbit`, `State`, `OrbitArray` and `StateArray` objects and allows some kind of propagation on all of them. 
# 
# This notebook uses type annotations and run-time type checks - for testing only. This stuff is not going into `poliastro` ...

# ## Import

# In[1]:


from datetime import datetime, timedelta
from math import atan2, sqrt
from random import random
from typing import Union

import astropy.units as u
from astropy.time import Time, TimeDelta
import numpy as np
from typeguard import typechecked


# ## Data
# 
# We need some test data, i.e. an array of time stamps and some initial points in 2D space.

# In[2]:


START = Time('2022-01-01 00:00', scale = 'ut1')
STOP = Time('2022-01-10 00:00', scale = 'ut1')
STEP_LENGTH = 1 << u.d

TIMES = Time(np.arange(START, STOP, STEP_LENGTH))
TIMES


# In[3]:


FACTOR = 10
POINTS_LENGTH = 5

POINTS = [
    {'x': random() * FACTOR << u.km, 'y': random() * FACTOR << u.km}
    for _ in range(POINTS_LENGTH)
]
POINTS


# ## The Baseline: "Scalar" `Orbit` and `State` Classes
# 
# Differing from `poliastro`, the `epoch` is moved into the `State` class (and its decedents). This could potentially allow some interesting optimizations and extensions when implementing a `StateArray` class because it centralizes all relevant data into one place.

# In[4]:


@typechecked
class State:
    
    @u.quantity_input(x = u.km, y = u.km)
    def __init__(self, epoch: Time, x: u.Quantity, y: u.Quantity):
        assert epoch.ndim == x.ndim == y.ndim == 0
        self._epoch = epoch
        self._x = x
        self._y = y
    
    def __repr__(self) -> str:
        return f'<State epoch={self._epoch} x={self._x} y={self._y}>'
    
    def __eq__(self, other) -> bool:
        return self._epoch == other.epoch and self._x == other.x and self._y == other.y
    
    @property
    def epoch(self) -> Time:
        return self._epoch
    
    @epoch.setter
    def epoch(self, value: Time):
        assert value.ndim == 0
        self._epoch = value
    
    @property
    def x(self) -> u.Quantity:
        return self._x
    
    @x.setter
    @u.quantity_input(value = u.km)
    def x(self, value: u.Quantity):
        assert value.ndim == 0
        self._x = value
    
    @property
    def y(self) -> u.Quantity:
        return self._y
    
    @y.setter
    @u.quantity_input(value = u.km)
    def y(self, value: u.Quantity):
        assert value.ndim == 0
        self._y = value
    
    def to_polar(self) -> tuple[u.Quantity, u.Quantity]:
        return (
            sqrt(self._x.to_value(u.km) ** 2 + self._y.to_value(u.km) ** 2) << u.km,
            atan2(self._y.to_value(u.km), self._x.to_value(u.km)) << u.rad,
        )


# We can set up a state and play with it:

# In[5]:


a = State(Time('2022-01-01 00:00', scale = 'ut1'), 2.0 << u.km, 3.0 << u.km)
print(a)
print(a.to_polar())


# Next up is the `Orbit` class. It is more or less a simple wrapper around the state from a data perspective. Otherwise it would keep virtually all of its established methods. For simplicity in this PoC, almost all methods have been stripped and it directly includes its "propagator" code. Notice that the `propagate` method can handle both time scalars and time arrays, the latter resulting in an `OrbitArray` object. Also notice that there is a new `inplace` argument, making the orbit object as well as its state object (more explicitly) mutable.

# In[6]:


@typechecked
class Orbit:

    def __init__(self, state: State):
        self._state = state

    def __repr__(self) -> str:
        return f'<Orbit epoch={self._state.epoch} x={self._state.x} y={self._state.y}>'

    def propagate(self, timedelta: TimeDelta, inplace: bool = False): # HACK ignore the check ... 'Union[Orbit, OrbitArray]'
        
        if timedelta.ndim == 0:
        
            days = timedelta.to_value(u.d)
            dx = (random() - 0.5) * days << u.km
            dy = (random() - 0.5) * days << u.km
            
            if inplace:
                self._state.epoch += timedelta
                self._state.x += dx
                self._state.y += dy
                return self
            else:
                return type(self)(State(
                    epoch = self._state.epoch + timedelta,
                    x = self._state.x + dx,
                    y = self._state.y + dy,
                ))
        
        else:
            
            assert not inplace

            days = timedelta.to_value(u.d)
            dx = ((np.random.random(timedelta.size) - 0.5) << u.km).reshape(timedelta.shape) * days
            dy = ((np.random.random(timedelta.size) - 0.5) << u.km).reshape(timedelta.shape) * days
            
            return OrbitArray(StateArray(
                epoch = self._state.epoch + timedelta,
                x = self._state.x + dx,
                y = self._state.y + dy,
            ))

    @property
    def state(self) -> State:
        return self._state


# We can set up an orbit (via a state) and play with it. In the first example, the `propagate` generates a new `Orbit` object as `poliastro` currently does:

# In[7]:


a = Orbit(State(Time('2022-01-01 00:00', scale = 'ut1'), 2.0 << u.km, 3.0 << u.km))
print(a)
b = a.propagate(TimeDelta(1 << u.d), inplace = False)
print(b)
print(a is b)
print(a == b)


# The second example does the conversion in place but returns the object itself for consistency and convenience:

# In[8]:


c = a.propagate(TimeDelta(1 << u.d), inplace = True)
print(c)
print(a is c)
print(a == c)


# ## New Territory: The `OrbitArray` and `StateArray` Classes
# 
# Time to look at the array types. It has a number of methods and properties that make it look like a `numpy.ndarray`. Like its `numpy` counterpart, it is in fact actually n-dimensional. In theory, longer-term, it could actually be an object exposing an [array interface](https://numpy.org/doc/stable/reference/arrays.interface.html). By contrast, it could also be something else entirely and the one place which is responsible for stripping units (and keeping track of them separately). Anyway, a simple first sketch that at least allows most common types of slicing and falls back to "scalar" `State` objects in the right places etc:

# In[9]:


@typechecked
class StateArray:

    @u.quantity_input(x = u.km, y = u.km)
    def __init__(self, epoch: Time, x: u.Quantity, y: u.Quantity):
        assert epoch.shape == x.shape == y.shape
        self._epoch = epoch
        self._x = x
        self._y = y

    def __repr__(self) -> str:
        return  (
            f'<StateArray shape={self.shape} value=[\n'
            + '\n'.join([
                f' (epoch={state.epoch} x={state.x} y={state.y}),'
                for state in self.reshape(self.size)
            ])
            + '\n]>'
        )
    
    def __getitem__(self, idx) -> 'Union[StateArray, State]':
        target = type(self)(
            epoch = self._epoch[idx],
            x = self._x[idx],
            y = self._y[idx],
        )
        if np.squeeze(target.epoch).ndim == 0 and (
            isinstance(idx, int) or (
                isinstance(idx, tuple) and all(isinstance(item, int) for item in idx)
            )
        ):
            return State(
                epoch = target.epoch,
                x = target.x,
                y = target.y,
            )
        return target

    def reshape(self, *args) -> 'StateArray':
        return type(self)(
            epoch = self._epoch.reshape(*args),
            x = self._x.reshape(*args),
            y = self._y.reshape(*args),
        )

    @property
    def epoch(self) -> Time:
        return self._epoch

    @property
    def x(self) -> u.Quantity:
        return self._x

    @property
    def y(self) -> u.Quantity:
        return self._y

    @property
    def ndim(self):
        return self._epoch.ndim

    @property
    def size(self):
        return self._epoch.size

    @property
    def shape(self):
        return self._epoch.shape

    def to_polar(self) -> tuple[u.Quantity, u.Quantity]:
        return np.sqrt(self._x ** 2 + self._y ** 2), np.arctan2(self._y, self._x)

    @classmethod
    def from_states(cls, states):
        return cls(
            epoch = Time([state.epoch for state in states]),
            x = u.Quantity([state.x for state in states], u.km),
            y = u.Quantity([state.y for state in states], u.km),
        )


# We can initialize a list of `State` objects and use them to create a `StateArray` object before playing with it:

# In[10]:


states = [
    State(time, **position)
    for time, position in zip(TIMES, POINTS)
]
states


# In[11]:


statearray = StateArray.from_states(states)
a = statearray[:4].reshape(2, 2)
a


# In[12]:


a[1, 1]


# In[13]:


statearray.to_polar()


# Just like the `Orbit` class is a wrapper around the `State` class in many ways, the `Orbit` array mostly wraps `StateArray`. Its `propagate` method could theoretically allow standard `numpy` broadcasting rules depending on shape of the `timedelta` parameter. In the following example, for simplicity, only scalar values for `timedelta` or values for `timedelta` with a shape matching the shape of the `OrbitArray` are allowed:

# In[14]:


@typechecked
class OrbitArray:
    
    def __init__(self, statearray: StateArray):
        self._statearray = statearray
    
    def __repr__(self) -> str:
        return  (
            f'<OrbitArray shape={self.shape} value=[\n'
            + '\n'.join([
                f' (epoch={state.epoch} x={state.x} y={state.y}),'
                for state in self._statearray.reshape(np.multiply.reduce(self._statearray.shape))
            ])
            + '\n]>'
        )
    
    def __getitem__(self, idx) -> 'Union[OrbitArray, Orbit]':
        target = self._statearray[idx]
        if isinstance(target, State):
            return Orbit(state = target)
        return OrbitArray(statearray = target)

    def propagate(self, timedelta: TimeDelta, inplace: bool = False) -> 'OrbitArray':
        
        if timedelta.ndim == 0:
            timedelta = np.repeat(timedelta.to_value(u.d), self.size).reshape(self.shape) << u.d
        else: # TODO allow better broadcasting logic
            assert timedelta.shape == self.shape
        
        days = timedelta.to_value(u.d)
        dx = ((np.random.random(self.size) - 0.5) << u.km).reshape(self.shape) * days
        dy = ((np.random.random(self.size) - 0.5) << u.km).reshape(self.shape) * days
        
        if inplace:
            self._statearray.epoch[:] += timedelta
            self._statearray.x[:] += dx
            self._statearray.y[:] += dy
            return self
        else:
            return type(self)(StateArray(
                epoch = self._statearray.epoch + timedelta,
                x = self._statearray.x + dx,
                y = self._statearray.y + dy,
            ))

    def reshape(self, *args) -> 'OrbitArray':
        return type(self)(statearray = self._statearray.reshape(*args))

    @property
    def statearray(self) -> StateArray:
        return self._statearray

    @property
    def ndim(self):
        return self._statearray.ndim

    @property
    def shape(self):
        return self._statearray.shape

    @property
    def size(self):
        return self._statearray.size

    @classmethod
    def from_orbits(cls, orbits):
        return cls(
            statearray = StateArray.from_states([orbit.state for orbit in orbits]),
        )


# Using our previously created `StateArray` object, we can now play with the `OrbitArray` class:

# In[15]:


a = OrbitArray(statearray)
a


# In[16]:


b = a[:4].reshape(2, 2)
b


# In[17]:


c = b.propagate(TimeDelta(7 << u.d))
b, c


# Because a single "scalar" `Orbit` object can be propagated to many points in time, we can now actually test it:

# In[18]:


a = Orbit(State(Time('2022-01-01 00:00', scale = 'ut1'), 2.0 << u.km, 3.0 << u.km))
a


# In[19]:


b = a.propagate(TimeDelta(np.arange(1, 7).reshape(2, 3) << u.d))
b


# ## Further thinking
# 
# - [Structured arrays](https://numpy.org/doc/stable/user/basics.rec.html) can actually be n-dimensional. This could be a foundation for `StateArray` classes, opening the back door for alternative libraries which are exposing an [array interface](https://numpy.org/doc/stable/reference/arrays.interface.html). In a bigger picture, the state array could be the place where stuff like `cupy.ndarray` or `dask.array` are transparently accepted.
# - [Dask and xarray vs Astropy](https://github.com/astropy/astropy/issues/12600) - this would not be the first time people are trying it. To enable it, one actually needs to take care of the Units somehow. Since we need to strip them anyway for `numba`, ideally in one central place, `StateArray` could be the one place where it happens.
# - Numpy `ufunc`s can apparently (transparently) be used on `cupy` arrays, experimental though, see this [issue](https://github.com/cupy/cupy/issues/2011). This is not a general rule, `cupy` needs to have its own version of the function. The beauty is that one does not need to write different code - `cupy` automatically patches `numpy` function calls of they are applied to a `cupy` array.
# - Dask arrays are known to work with `numba` via `gufunc`s, see [here](https://docs.dask.org/en/latest/generated/dask.array.gufunc.apply_gufunc.html) and [here](https://examples.dask.org/applications/stencils-with-numba.html).
# - Just a note: `astropy`'s [Time](https://docs.astropy.org/en/stable/time/index.html) arrays are a little inconsistent or flexible, depending on the point of view. They can reference time stamps as strings as well as `datetime` objects. Not sure if this is an issue here.

# In[ ]:




