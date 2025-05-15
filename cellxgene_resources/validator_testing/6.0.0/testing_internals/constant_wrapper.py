"""
Wrapper class to use dict[str[str]] to turn the str values into an iterable/list
that agrees with pytest.

TEST_CONST = TestConstant({
    "string_key_descriptor": "string test value",
    "string_key_descriptor2": "string test value2",
})

Can use standard[key] = value on instance to get back value wrapped in list
.values() method will return list of original input values

Can add more later if it is useful to have different values beyond strings
"""

class TestConstant:
    def __init__(self, test_dict):
        self.test_dict = test_dict
        self._test_data_dict = None
        self._values_to_list()

    def _values_to_list(self):
        assert all(isinstance(key, str) for key in self.test_dict), "Input keys must be strings"
        values_dict = {key: [value] for key, value in self.test_dict.items()}
        assert all(isinstance(value, list) for value in values_dict.values())
        self._test_data_dict = values_dict

    def __getitem__(self, key):
        return self._test_data_dict[key]

    def values(self):
        return self.test_dict.values()
