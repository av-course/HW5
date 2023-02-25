"""

Instructions:

Answer the following questions on paper (or typeset your solutions using LaTeX). If handwriting solutions, everything must be very neat and legible.

************************************************
Problem 1

Consider a more detailed vehicle model than we discussed in class. This model is known as the dynamic bicycle model, and is a more accurate model of how vehicles move. 

xₖ = [p₁
      p₂
      θ
      v₁
      v₂
      ω]

uₖ = [a
      δ]


xₖ₊₁ = xₖ + Δ * [cos(θ)*v₁ - sin(θ)*v₂
                 sin(θ)*v₁ + cos(θ)*v₂
                 ω
                 a + v₂*ω
                 -v₁ * ω + 2 v₂ / v₁ - δ
                 2 * ω / v₁ - δ]

Consider a function F(xₖ, uₖ), which produces xₖ₊₁

a) Derive the Jacobian of F w.r.t. xₖ. This should be a 6x6 matrix, which potentially involves terms from xₖ, uₖ, and Δ.
b) Derive the Jacobian of F w.r.t. uₖ. This should be a 6x2 matrix, which potentially involves terms from xₖ, uₖ, and Δ.

************************************************
Problem 2

Consider the problem of approximating a centerline (which runs between two road boundaries, represented by known points), as an implicit quadratic curve.

Assume we have three points along each road boundary, lᵢ, 1≤i≤3, and rᵢ, 1≤i≤3.

               *l₃
        *l₂
    *l₁  
                   *r₃
        *r₂
    *r₁

We wish to encode the centerline as the implicit curve satisfying f(p) = 0, where f(p) := a * (p₁² + p₂²) + b₁*p₁ + b₂*p₂.

To find the parameters a, b₁, and b₂ which best approximate the centerline, we can solve the following optimization problem:

min        Σᵢ ( (f(rᵢ) - 1)² + (f(lᵢ) + 1)² )
a,b₁,b₂


Express the above optimization as a standard quadratic program, i.e. a program of the form

min                      [a
a,b₁,b₂   [a b₁ b₂]* Q *  b₁   + [a b₁ b₂] * q
                          b₂]

for an appropriate matrix Q and vector q, which you must find.

************************************************
Problem 3

Assume p(xₖ₋₁|z₁, ..., zₖ₋₁), p(xₖ|xₖ₋₁), and p(zₖ|xₖ) are known.

Show how p(xₖ|z₁, ..., zₖ) can be expressed in terms of the known distributions.


