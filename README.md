# Scripts and functions for interacting with VLASS data

## Installation

1. Clone repo and `cd` into repo directory.

2. `python setup.py install`

## Use

In Python:

```
> from vlass_tools import get_vlass
> co = (<your_ra>, <your_dec>)
> get_vlass.get_tilename(co)
```
