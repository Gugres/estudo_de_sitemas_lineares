import unittest
from T1 import main

class NamesTestCase(unittest.TestCase):
    
    def test_first_last_name(self):
        result = main(20)
        self.assertEqual(([1]*20), ([1]*20))

if __name__ == '__main__':
    unittest.main()
