function result = isPositiveInteger(value)
    result = isnumeric(value) && isreal(value) && value > 0 && mod(value, 1) == 0;
end