"""CSV export must not undo the upstream 200 Hz clock correction."""
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class CsvTimestampTests(unittest.TestCase):
    def test_full_duration_clock_round_trips_and_preserves_column_format(self):
        source = r'''
#include "util/CsvTimestamp.h"
#include <sstream>
#include <string>
int main() {
    for (int i = 0; i < 240000; ++i) {
        const double time = double(i) / 200.0;
        std::ostringstream out;
        out << std::fixed << std::setprecision(2);
        ocean_imu::sim::write_csv_time(out, time) << ',' << 1.234567;
        const auto text = out.str();
        const auto comma = text.find(',');
        if (std::stod(text.substr(0, comma)) != time) return 1;
        if (text.substr(comma + 1) != "1.23") return 2;
    }
}
'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / 'test.cpp').write_text(source)
            subprocess.run(['g++', '-std=c++17', '-O2', '-I', str(ROOT / 'src'),
                            str(path / 'test.cpp'), '-o', str(path / 'test')], check=True,
                           capture_output=True, text=True)
            subprocess.run([str(path / 'test')], check=True, capture_output=True, text=True)


if __name__ == '__main__':
    unittest.main()
