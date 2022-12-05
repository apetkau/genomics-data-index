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

import logging

logger = logging.getLogger(__name__)

try:
    from ete3 import TreeStyle, NodeStyle, Face, RectFace, CircleFace, TextFace
except ImportError as e:
    logger.warning("Could not import ete3 package. Visualization of dendrograms is unavailable.", e)

    import os

    create_mock_classes = True

    # Try to set 'QT_QPA_PLATFORM' as specified in https://github.com/etetoolkit/ete/issues/500
    if 'QT_QPA_PLATFORM' not in os.environ:
        logger.warning("QT_QPA_PLATFORM unset. Attempting to set QT_QPA_PLATFORM='offscreen' and import ete3 package")
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'

        try:
            from ete3 import TreeStyle, NodeStyle, Face, RectFace, CircleFace, TextFace
            logger.warning("Could not import ete3 package after adjusting QT_QPA_PLATFORM. Visualization of dendrograms is unavailable")
            create_mock_classes = False
        except ImportError:
            pass

    # If cannot import appropriate modules, create dummy/mock objects
    if create_mock_classes:
        error_msg = ("Could not properly import {}. Error message: [{}]. "
                "If the ete3 package is found, then this is likely due to a missing or improperly installed X server, "
                "which is required for graphical functionality "
                "within the ETE toolkit (see <https://github.com/etetoolkit/ete/issues/101>). "
                "Please either install an X server or attempt to run the application within a virtual framebuffer (like 'xvfb')")

        class TreeStyle:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try:
                    from ete3 import TreeStyle
                except ImportError as e:
                    raise Exception(error_msg.format('TreeStyle', str(e)), e)

                raise Exception('This should not occur')


        class NodeStyle:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try:
                    from ete3 import NodeStyle
                except ImportError as e:
                    raise Exception(error_msg.format('NodeStyle', str(e)), e)

                raise Exception('This should not occur')


        class Face:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try:
                    from ete3 import Face
                except ImportError as e:
                    raise Exception(error_msg.format('Face', str(e)), e)

                raise Exception('This should not occur')


        class RectFace:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try:
                    from ete3 import RectFace
                except ImportError as e:
                    raise Exception(error_msg.format('RectFace', str(e)), e)

                raise Exception('This should not occur')


        class CircleFace:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try:
                    from ete3 import CircleFace
                except ImportError as e:
                    raise Exception(error_msg.format('CircleFace', str(e)), e)

                raise Exception('This should not occur')


        class TextFace:

            # Force class to raise an exception on creation of an instance
            # so that errors are only returned to a user when this is used
            def __init__(cls):
                try:
                    from ete3 import TextFace
                except ImportError as e:
                    raise Exception(error_msg.format('TextFace', str(e)), e)

                raise Exception('This should not occur')
