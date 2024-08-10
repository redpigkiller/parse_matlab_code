import unittest
from pathlib import Path

from parse_matlab_code.core.tokenizer import Tokenizer
from parse_matlab_code.core.parser_error import *


class TestTokenizerr(unittest.TestCase):

    def setUp(self) -> None:
        self.sample_code_pathes = [path for path in Path('./test_data').rglob('*.m')]

    def test_tokenizer(self):
        tokenizer = Tokenizer()

        for index, sample_code_path in enumerate(self.sample_code_pathes):
            with open(sample_code_path, mode='r', encoding='utf-8') as f:
                code = f.read()
            print(f"Running tokenizer test {index + 1}/{len(self.sample_code_pathes)}...")

            with self.subTest(i=index):
                try:
                    tokens = tokenizer.tokenize(code)
                except Exception as e:
                    self.fail(f"During tokenizing {sample_code_path}\nTokenizer raised an exception: {e}")

if __name__ == '__main__':
    unittest.main()
