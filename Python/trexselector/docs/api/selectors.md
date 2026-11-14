# Selectors and controls

## Selectors

The seven selectors of the package root, the classes that the guide pages use.

```{eval-rst}
.. autoclass:: trexselector.TRexSelector

.. autoclass:: trexselector.TRexDASelector

.. autoclass:: trexselector.TRexGVSSelector

.. autoclass:: trexselector.TRexTikhonovSelector

.. autoclass:: trexselector.TRexScreeningSelector

.. autoclass:: trexselector.TRexBiobankScreeningSelector

.. autoclass:: trexselector.TRexSPCASelector
```

## Controls, results and core classes

The control and result classes and the enums of the selectors, as re-exported at the package root. The
module also holds the core class of every selector, the binding of the C++ class with its getter
methods, which the selector of the package root wraps.

```{eval-rst}
.. automodule:: trexselector.trex_selector_methods.ad
   :imported-members:
```
