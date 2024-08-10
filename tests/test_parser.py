import unittest
from pathlib import Path

from parse_matlab_code.core.parser import Parser
from parse_matlab_code.core.parser_error import *


class TestTokenizerr(unittest.TestCase):

    def setUp(self) -> None:
        self.sample_code_pathes = [path for path in Path('./test_data').rglob('*.m')]

    def test_tokenizer(self):
        parser = Parser()

        for index, sample_code_path in enumerate(self.sample_code_pathes):
            with open(sample_code_path, mode='r', encoding='utf-8') as f:
                code = f.read()
            print(f"Running parser test {index + 1}/{len(self.sample_code_pathes)}...")
            
            with self.subTest(i=index):
                try:
                    parse_tree = parser.parse(code)
                except Exception as e:
                    self.fail(f"During parsing {sample_code_path}\nParser raised an exception: {e}")

if __name__ == '__main__':
    unittest.main()
