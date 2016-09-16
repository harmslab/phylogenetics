Handlers
========

Most of phylogenetics is gathering sequence data, converting it to different file
formats, tracking the output at each step along the process, and trying to remember
what you do six months ago. All of this is very time consuming. The `phylogenetics`
python package takes care of this internally using objects known as `Handlers`.
They are are simple, useful templates to store sequence data and track all history
of that object. Everytime a change is made, a `Handler` stores it in its history.
Further, we've designed a `HandlerContainer` to manage these objects and track
whole sets through a phylogenetics project. All the pieces you see in the `phylogenetics`
package are subclasses of these two objects.

If the ``phylogenetics`` package does not have a ``Handler`` that you need, you
can easily subclass and write your own custom one!

Basic Examples
--------------

The example below is a basic custom ``Handler`` subclass. This allows you to add
methods that are specific to the handler you are creating.

.. code::

    from phylogenetics.handlers import Handler

    class CustomHandler(Handler):
        """Put your docstring here.
        """
        def __init__(self, **kwargs):
            super(CustomHandler, self).__init__(**kwargs)

The example below is a custom ``HandlerContainer`` subclass. This can container multiple
Handlers of a specific type. You can add your methods for manipulating whole
groups of handlers as well.

.. code::

    from phylogenetics.handlers import HandlerContainer

    class CustomContainer(HandlerContainer):
        """Put your docstring here.
        """
        def __init__(self, **kwargs):
            super(CustomHandler, self).__init__(**kwargs)

        @property
        def _prefix(self):
            return "ITEM"

        @property
        def _child_type(self)
            """The Handler object subclass passed into this Container. Must be
            a contained in a list. Can take in many types of Handlers
            """
            return [CustomHandler]
