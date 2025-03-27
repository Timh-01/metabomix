# Goes in add.py

def add(number1: int | float, number2: int | float) -> int:
    """
    Integer addition, if floats are provided they will be first
    converted to integers by rounding down

    Examples:
        >>> add(1, 2)
        3

        >>> add(2.3, 4.5)
        6

    Args:
        number1 (int | float): first number for addition
        number2 (int | float): second number for addition

    Returns:
        int: sum of (possibly rounded down) inputs
    """
    return int(number1) + int(number2)