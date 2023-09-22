# Integer Factorization via unit groups

This program implements a kind of algebraic group factorization, using two families of groups.


**The group of solutions to $x^2 + dy^2 = 1$**

The group operation is defined like multiplication in the algebraic number field $Q(\sqrt{-d})$. These groups have order $p \pm 1$, depending on the choice of $d$ and whether $-1$ is a square mod p.

**The group of solutions to $x^2 + d(y^2 + z^2 + w^2) = 1$**

Here, the group operation is similar to quaternion multiplication.

$$
(x_1, y_1, z_1, w_1) * (x_2, y_2, z_2, w_2) =
$$

$$
(x_1 x_2 - d(y_1 y_2 + z_1 z_2 + w_1 w_2), \space x_1 y_2 + y_1 x_2 + z_1 w_2 - w_1 z_2,
$$

$$
x_1 z_2 + z_1 x_2 + w_1 y_2 - y_1 w_2, \space x_1 w_2 + w_1 x_2 + y_1 z_2 - z_1 y_2)
$$

Possible orders here are $p^3 \pm p$.
