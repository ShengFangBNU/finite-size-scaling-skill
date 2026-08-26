"""sys.path bootstrap so scripts run as ``python scripts/<name>.py`` can
import the ``fss`` package from the repository root."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
