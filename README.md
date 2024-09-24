# bigInt

This library implements main operations with long integer numbers. Inside of bigInt.h file you will find a bigInt class and some more functions. Inside of the class, the numbers are stored base 2^bits (where bits = 10) leat significant digit first. This means that number X will need O(log X / log base) = O(log X / 10) memory (let`s call n = log X / 10 -- number of digits of a number). bigInt:digits is array of digits of a number, bigInt:sign is 1 iff number is negative.

Here are operations implemented with numbers:

## Сonversion from int/long long

Simply convert given number to needed base.

Complexity: O(n).

## Inversion(unary -)

Just swaps the sign bit.

Complexity: O(1).

## Сomparison(<, >, ==, <=, >= )

Going through numbers digit by digit from the most significant ones, we find the first different digit.

Complexity: O(n).

## Addition(+, +=, ++)

The classic digit-by-digit addition with carrying in given base. 

Complexity: O(n).

## Subtraction(-, -=, --)

Similar to addition algorithm that works in O(n).

## Bitwise operations(&, &=, |, |=, ^, ^=)

Because we store numbers in base which is power of two, we can just do operations digit by digit and it will give us the right result.
 
Complexity: O(n)

## Bit shifts (<<, <<=, >>, >>=)

If we shift number by k bits to left, first we delete first (k / bits) digits, and after that we need to move (k % bits) from each digit to the previous one.

Similar, with right shift, we add (k / bits) zeros to the beginning and then move last (k % bits) last bits of each number to the next one.

Complexity: O(n).

## Multiplication (*, *=)

Multiplication is tricky because the usual method runs in O(n ^ 2), which is too long. I also tried to implement the karatsuba algorithm, which finds the product of two numbers using the divide-and-conquer method in O(n ^ (log_2(3))) time. This is much better than quadratic multiplication, but in practice it also takes quite a long time. So I used a faster multiplication algorithm with Fast Fourier Transform(FFT).

The key idea is to create polynomials whose coefficients will be the digits of the number, multiply these two polynomials and reduce the product back to the form of a number. The non-trivial part is to multiply two polynomials. In order to do this, let's redefine the polynomials in a different form -- the values ​​of the polynomials at a sufficient number of some points. In fact, this is what the FFT algorithm can do. It finds the value of a polynomial at the points w^1, w^2, ... . Where w is a primitive root of unity. After that we just multiply values of polynomials in the points and get the resulting polynomial. Now we can find transform a polynomial to standart form using inverse Fourier Transform.

The complexity of FFT is O(n * log n), all the other transformations is done in O(n).

Complexity: O(n * log n)

## Integer division (/, /=)

To begin with, we need to consider cases of division of negative numbers. Consider the case of a positive divisor and a dividend. In order to find the ratio a / b, first find the reciprocal d = 1 / b, and then the answer will be a * d. To find d, let's use Newton's method. Newton's method is used for an approximate search for the zeros of a function f(x). In our case, this corresponds to finding the zero of the function f(x) = b * x - 1. Then the transition to approximate the value of d is x1 = x * (2 - x * d). With the correct choice of the value of x0, each such transition will improve the value of d by a factor of two, so we will need to perform O(log n) operations. 

The problem with this method is non-integer d and intermediate results. This can be solved by finding d' = 2^m * d for some m (let's use m = popcount(a) + popcount(b)). Then the transition will look like x' -> (2 * x' * (2^m) - x' * d) / (2 ^ m) = (2 * x' * (2^m) - x' * d) >> m. We see that operations do not use non-integer operations, therefore such an algorithm is suitable for us. So we can find answer as (a * x') / (2^m). As good starting value is x' = a - b.

So we will do the O(log n) steps, each step takes O(n * log n).

Complexity: O(n * log^2 n).

## Mod (%, %=)

Just used the division algorithm and return a % b = a - (a / b) * b.

Complexity: O(n * log^2 n).

## Power (pow)

Classic binary exponentation, where we each time divide the exponent by two.

Complexity: O(n log^2 n) - O(log n) multiplications.

## Output in another base (change_base, to_bin, to_hex, to_decimal)

change_base returns array of digits of the number in given base, most significant digit first.
to_bin, to_decimal and to_hex return strings which is the number in binary, decimal and hexadecimal representation respectively.

Let's see how change_base works when changing base to _base. It is done with divide-and-conquer technique. We recursivly find first half and second halfs of digits base _base, concatenate them and return the answer. More precisely, for number x, m is equal to number of digits of x base _base (WLOG n is a power of 2). We find the result recursively. Let`s say x = L + (_base ^ (m / 2)) * R. So we just need to concatenate L.change_base() and R.change_base() which would be a result.

Function change_base uses additional function change_base_rec, which is private.

On each level of recursion the m is divided by 2, so we need a O(log n) steps to find a answer.

Complexity: O(n log^3 n).

## Input from another base (convert, initialisation from a decimal string)

The function convert takes an array of digits d in base _base (least digits first) and returns a bigInt number, which is equal to given.

It is also done with divide-and-conquer. WLOG |d| = power of two. We recursively find L = convert(d[:|d| / 2]) and R = convert(d[|d| / 2:]), so the answer can be found as R + L * (_base * (|d| / 2)).

Function convert uses a additional function convert_rec.

The time complexity analisis is analogical to change_base.

Complexity: O(n log^3 n).

Some other functions:

## pop_count

Returns the number of bits of a number.

Complexity: O(n).

## gcd

The standart Euclidean algorithm.

COmplexity: O(n^2 * log n).



