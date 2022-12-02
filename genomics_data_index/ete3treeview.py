# This file encapsulates the imports from the ETE Toolkit which rely on graphics.
# Importing these packages from the ETE Toolkit requires an X server to be installed
# (see https://github.com/etetoolkit/ete/issues/500).
# I wished to allow the Genomics Data Index (GDI) to be installed and used on machines without 
# an X server, and hence there is a complicated bit of importing of Python packages here
# which will replace the corresponding graphics classes (e.g., TreeStyle) from the ETE Toolkit
# with a dummy class that raises an exception when it is first used.
# A better approach for the future here would be to separate out any graphics functionality 
# into a separate Python package which can be skipped for those who do not wish to use 
# GDI in a graphical environment.

try:
    from ete3 import TreeStyle, NodeStyle, Face, RectFace, CircleFace, TextFace
except ImportError:
    import os

    create_mock_classes = True

    # Try to set 'QT_QPA_PLATFORM' as specified in https://github.com/etetoolkit/ete/issues/500
    if 'QT_QPA_PLATFORM' not in os.environ:
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'

        try:
            from ete3 import TreeStyle, NodeStyle, Face, RectFace, CircleFace, TextFace
            create_mock_classes = False
        except ImportError:
            pass

    # If cannot import appropriate modules, create dummy/mock objects
    if create_mock_classes:
        error_msg = ("Could not properly import {}. Error message {}.\n"
                "This is likely due to a missing or improperly installed X server, which is required for graphical functionality "
                "within the ETE toolkit (see <https://github.com/etetoolkit/ete/issues/101>).\n"
                "Please either install an X server or attempt to run the application within a virtual framebuffer (like 'xvfb')")

        class TreeStyle:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try
                    from ete3 import TreeStyle as TreeStyle1
                except ImportError as e:
                    raise Exception(error_msg.format('TreeStyle', str(e)), e)

                raise Exception('This should not occur')


        class NodeStyle:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try
                    from ete3 import NodeStyle as NodeStyle1
                except ImportError as e:
                    raise Exception(error_msg.format('NodeStyle', str(e)), e)

                raise Exception('This should not occur')


        class Face:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try
                    from ete3 import Face as Face1
                except ImportError as e:
                    raise Exception(error_msg.format('Face', str(e)), e)

                raise Exception('This should not occur')


        class RectFace:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try
                    from ete3 import RectFace as RectFace1
                except ImportError as e:
                    raise Exception(error_msg.format('RectFace', str(e)), e)

                raise Exception('This should not occur')


        class CircleFace:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try
                    from ete3 import CircleFace as CircleFace1
                except ImportError as e:
                    raise Exception(error_msg.format('CircleFace', str(e)), e)

                raise Exception('This should not occur')


        class TextFace:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try
                    from ete3 import TextFace as TextFace1
                except ImportError as e:
                    raise Exception(error_msg.format('TextFace', str(e)), e)

                raise Exception('This should not occur')
