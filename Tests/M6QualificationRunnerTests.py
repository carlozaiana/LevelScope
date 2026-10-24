"""Fail closed on incomplete/parity-violating M6 capture evidence."""
import pathlib
import struct
import tempfile
import unittest
from compare_m6_reports import compare_capture, compare


class CaptureTests(unittest.TestCase):
    def test_identity(self):
        data = struct.pack('<4f', 0, -.3, 1.25, -2)
        self.assertEqual(compare_capture(data, data, True), 0)

    def test_one_ulp_nonzero_only(self):
        a, b = struct.pack('<I', 0x3f800000), struct.pack('<I', 0x3f800001)
        self.assertEqual(compare_capture(a, b, False), 1)
        with self.assertRaises(ValueError): compare_capture(a, b, True)
        with self.assertRaises(ValueError): compare_capture(a, struct.pack('<I', 0x3f800002), False)

    def test_negative_and_nonfinite(self):
        self.assertEqual(compare_capture(struct.pack('<I', 0xbf800000), struct.pack('<I', 0xbf800001), False), 1)
        with self.assertRaises(ValueError): compare_capture(struct.pack('<f', float('nan')), struct.pack('<f', 0), False)

    def test_missing_and_truncated(self):
        with self.assertRaises(ValueError): compare_capture(b'\0', b'\0', False)
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaises(OSError): compare(pathlib.Path(root), pathlib.Path(root))


if __name__ == '__main__':
    unittest.main()
