import sympy as sp
from sympy.calculus.util import continuous_domain
import numpy as np
import matplotlib.pyplot as mpl

x = sp.symbols("x")

class PPF:
    non_linear_list = []
    concavity_list = []
    max_power_list = []

    def __init__(self, function, number, x_max, first_derivative,
                 second_derivative, is_concave, rounding, max_power):
        self.function = function
        self.number = number
        self.oc = first_derivative * -1
        self.x_max = x_max
        self.y_max = function.subs(x, 0)
        self.is_concave = is_concave
        self.rounding = rounding
        self.max_power = max_power

        if sp.simplify(second_derivative) != 0:
            PPF.non_linear_list.append(self)

        PPF.concavity_list.append(self.is_concave)
        PPF.max_power_list.append(self.max_power)


def EvaluateRounding(function):
    monomials = function.as_ordered_terms()
    powers = []

    for monomial in monomials:
        monomial = str(monomial)
        positions = []

        if monomial.find("x") != -1:
            for position in range(len(monomial)):
                if monomial[position] == "*":
                    if monomial[position + 1] == "*":
                        positions.append(position)

            power = float(1)

            for position in positions:
                string = ""
                current = position + 2
                while current != len(monomial) and monomial[current] != " ":
                    string += monomial[current]
                    current += 1
                power *= float(string)
            powers.append(power)

    count_of_non_integers = 0
    all_int = True
    non_integer_powers = []
    for power in powers:
        if power.is_integer() and count_of_non_integers == 0:
            all_int = True

        if not power.is_integer():
            all_int = False
            count_of_non_integers += 1
            non_integer_powers.append(power)

    rounding = True
    suitable = True
    max_power = max(powers)

    if all_int:
        if max_power < 3:
            rounding = False

    else:
        for power in non_integer_powers:
            if not (power / 0.5).is_integer():
                suitable = False
    return rounding, suitable, max_power


def RequestAndCheckPPF(num_of_functions_entered):
    try:
        function = sp.sympify(input("Enter the first PPF: ")) if num_of_functions_entered == 0\
            else sp.sympify(input("Enter the second PPF: "))
        for i in range(len(str(function))):
            if not (str(function)[i].isdigit() or str(function)[i] in ["+", "-", "/", "*", "x", " ", "."]):
                return RequestAndCheckPPF(num_of_functions_entered)
    except:
        return RequestAndCheckPPF(num_of_functions_entered)

    rounding, suitable, max_power = EvaluateRounding(function)
    if not suitable:
        print("The function include a non-integer power which does not end with .5")
        return RequestAndCheckPPF(num_of_functions_entered)

    # To check that a function has a positive solution
    function_equals_zero = sp.Eq(function, 0)
    solutions_1 = [root for root in sp.solve(function_equals_zero, x) if (root.is_real and root >= 0)]

    if not solutions_1:
        print(f"The function y = {function} does not have any positive solutions")
        return RequestAndCheckPPF(num_of_functions_entered)
    if rounding:
        x_max = solutions_1[0].evalf(3)
    else:
        x_max = float(solutions_1[0])
    non_rounded_x_max = solutions_1[0]
    # To check that a function is continuous on [0, x_max]:
    domain = continuous_domain(function, x, sp.Reals)
    interval = sp.Interval(0, x_max)
    check_whether_continuous = interval.issubset(domain)
    if not check_whether_continuous:
        print(f"The function y = {function} is not strictly defined on [0, x_max]")
        return RequestAndCheckPPF(num_of_functions_entered)

    # To check that a function has a continuous derivative sign on [0, x_max]:
    first_derivative = sp.diff(function, (x, 1))
    first_derivative_equals_zero = sp.Eq(first_derivative, 0)
    solutions_2 = sp.solve(first_derivative_equals_zero, x)
    solutions_2 = [root for root in solutions_2 if (root.is_real and root > 0 and root < non_rounded_x_max)]
    if solutions_2:
        print(f"The function y = {function} does not strictly decreases on [0, x_max]")
        return RequestAndCheckPPF(num_of_functions_entered)

    # To check that a function is decreasing on [0, x_max]:
    midpoint = x_max / 2
    first_derivative_at_midpoint = first_derivative.subs(x, midpoint)
    if first_derivative_at_midpoint >= 0:
        print(f"The function y = {function} does not strictly decreases on [0, x_max]")
        return RequestAndCheckPPF(num_of_functions_entered)

    # To check that a function is positive on [0, x_max]:
    # (From the conditions above)

    # To check for continuous concavity on [0, x_max]:
    second_derivative = sp.diff(function, (x, 2))
    second_derivative_equals_zero = sp.Eq(second_derivative, 0)
    solutions_3 = sp.solve(second_derivative_equals_zero, x)
    solutions_3 = [root for root in solutions_3 if (root.is_real and root > 0 and root < non_rounded_x_max)]

    if solutions_3:
        print(f"The PPF {function} does not have constant concavity on [0, x_max]")
        return RequestAndCheckPPF(num_of_functions_entered)

    # To check whether a function is concave or convex
    second_derivative_at_midpoint = second_derivative.subs(x, midpoint)
    if second_derivative_at_midpoint < 0:
        is_concave = True
    else:
        is_concave = False
    return function, x_max, first_derivative, second_derivative, is_concave, rounding, max_power


def Main():
    x = sp.sympify("x")

    rules_request = input("If you want to see the rules for a function to be a PPF, type Yes"
                          " Otherwise, anything else: ")

    if rules_request == "Yes":
        print("""
        0. Function should have the powers to be positive integers or ended by .5
        Example: y = 1 x4 suit; y = 1 x^0.3 does.

        A function should:
        1. Have at least one positive root (so that [0,x_max] can be determined).
        Example: y = 1/x does not suit.   y = 1-x does suit.

        For the interval [0,x_max] a function should:

        2. Be continuous.
        Example:
        y = 1-x does suit
        y = 1/x does not suit - not defined for x = 0

        3. Be non-negative.
        Example:
        y = 5-x does suit.
        y = x - 5 does not suit.

        4. Be non-increasing. In other words, having non-negative derivative sign (but since the function is decreasing,
        it can be equal to one only at 0 or x_max).
        Example:
        y = 5-x^2 does suit.
        y = sqrt(x) does not suit.

        5. Have constant concavity.
        Example:
        y = 1 - sqrt(x) does suit.
        y = 1-(x-1)^3 does not suit.
        """)

    print("A function should be entered as a function of x, without y = or f(x) before. for example: 1-x")

    function, x_max, first_derivative, second_derivative, is_concave, rounding, max_power = RequestAndCheckPPF(0)
    ppf_1 = PPF(function, 1, x_max, first_derivative, second_derivative, is_concave, rounding, max_power)
    function, x_max, first_derivative, second_derivative, is_concave, rounding, max_power = RequestAndCheckPPF(1)
    ppf_2 = PPF(function, 2, x_max, first_derivative, second_derivative, is_concave, rounding, max_power)
    concavity_count = 0

    for concavity in PPF.concavity_list:
        if concavity:
            concavity_count += 1
    try:
        if len(PPF.non_linear_list) == 2:
            # To determine the concavity configuration
            if concavity_count == 2:
                if max(PPF.max_power_list) <= 3:
                    NonLinear_1(ppf_1, ppf_2, x)
                else:
                    print("TThe combination is too complex to compute for now")
            elif concavity_count == 1:
                NonLinear_2(ppf_1, ppf_2, x)
            else:
                NonLinear_3(ppf_1, ppf_2, x)

        elif len(PPF.non_linear_list) == 1:
            # To determine which function is linear and which one is non-linear
            if PPF.non_linear_list[0].number == 1:
                non_linear = ppf_1
                linear = ppf_2
            else:
                non_linear = ppf_2
                linear = ppf_1
            # To determine the concavity configuration
            if concavity_count == 1:
                LinearNonLinear_1(non_linear, linear, x)
            else:
                LinearNonLinear_2(non_linear, linear, x)
        else:
            Linear(ppf_1, ppf_2, x)
    except:
        print("The combination is too complex to compute for now")


def UserInterface(output_list, function_list, boundary_list):
    x = sp.symbols("x")
    for bit in output_list:
        print(bit)
    print()
    boundary_list = sorted(list(set(boundary_list)))
    if boundary_list[0] == 0:
        boundary_list.pop(0)

    choice = input("Type Y if you want the graph to be plotted: ")
    if choice == "Y":
        boundary_list = sorted(list(set(boundary_list)))
        boundary_list = boundary_list if boundary_list[0] == 0 else [0] + boundary_list
        x_values = np.linspace(0, max(boundary_list), 1000)
        graph_boundaries = [(x_values >= boundary_list[i]) & (x_values < boundary_list[i + 1]) for i in
                            range(len(function_list))]
        function = []
        for bit in function_list:
            function.append(sp.lambdify(x, bit, "numpy"))
        y_vals = np.piecewise(x_values, graph_boundaries, function)
        mpl.plot(x_values, y_vals)
        mpl.xlabel("x")
        mpl.ylabel("y")
        mpl.title("Overall PPF")
        mpl.grid(True)
        mpl.show()

# Linear + Non-linear, non-linear is concave
def LinearNonLinear_1(non_linear, linear, x):
    output_list = []
    function_list = []
    boundary_list = []
    x_max = non_linear.x_max + linear.x_max
    x_max = x_max if type(x_max) == sp.Mul else round(float(x_max), 3)
    oc_non_linear_max = non_linear.oc.subs(x, non_linear.x_max)
    oc_non_linear_min = non_linear.oc.subs(x, 0)

    if linear.oc <= oc_non_linear_max and linear.oc >= oc_non_linear_min:
        equation = sp.Eq(linear.oc, non_linear.oc)
        equation_solution = sp.solve(equation, x)
        equation_solution = [solution for solution in equation_solution if (solution.is_real and solution >= 0)]
        boundary_point_1 = equation_solution[0]
        boundary_point_2 = boundary_point_1 + linear.x_max

        piece_1 = non_linear.function + linear.function.subs(x, 0)
        piece_2 = linear.function + piece_1.subs(x, boundary_point_1) - linear.y_max
        piece_2 = piece_2.subs(x, x - boundary_point_1)
        piece_3 = non_linear.function.subs(x, x - linear.x_max)
        boundary_point_1 = boundary_point_1 if type(boundary_point_1) == sp.Mul else round(float(boundary_point_1), 3)
        boundary_point_2 = boundary_point_1 if type(boundary_point_2) == sp.Mul else round(float(boundary_point_2), 3)
        output_list.append(f"y = {piece_1}, 0 <= x <= {boundary_point_1}")
        output_list.append(f"y = {piece_2}, {boundary_point_1} < x <= {boundary_point_2}")
        output_list.append(f"y = {piece_3}, {boundary_point_2} < x <= {x_max}")
        boundary_list.append(boundary_point_1), boundary_list.append(boundary_point_2), boundary_list.append(x_max)
        function_list.append(piece_1), function_list.append(piece_2), function_list.append(piece_3)


    elif linear.oc <= oc_non_linear_min:
        piece_1 = linear.function + non_linear.function.subs(x, 0)
        piece_2 = non_linear.function.subs(x, x - linear.x_max)
        output_list.append(f"y = {piece_1}, 0 <= x <= {linear.x_max}")
        output_list.append(f"y = {piece_2}, {linear.x_max} <= x <= {x_max}")
        boundary_list.append(linear.x_max), boundary_list.append(x_max)
        function_list.append(piece_1), function_list.append(piece_2)

    else:
        piece_1 = non_linear.function + linear.function.subs(x, 0)
        piece_2 = linear.function.subs(x, x - non_linear.x_max)
        non_linear.x_max = non_linear.x_max if type(non_linear.x_max) == sp.Mul else round(float(non_linear.x_max), 3)
        output_list.append(f"y = {piece_1}, 0 <= x <= {non_linear.x_max}")
        output_list.append(f"y = {piece_2}, {non_linear.x_max} <= x <= {x_max}")
        boundary_list.append(non_linear.x_max), boundary_list.append(x_max)
        function_list.append(piece_1), function_list.append(piece_2)

    return UserInterface(output_list, function_list, boundary_list)


# Two non-linears, both concave
def NonLinear_1(ppf_1, ppf_2, x):
    function_list = []
    output_list = []
    boundary_list = []
    oc_equality_x = sp.Eq(ppf_1.oc, ppf_2.oc)
    oc_equality_x_solution = sp.solve(oc_equality_x, x)

    if not oc_equality_x_solution:
        output_list.append(f"y = {ppf_1.function.subs(x, 0) + ppf_2.function}, 0 <= x <= {ppf_1.x_max + ppf_2.x_max}")
        function_list.append(ppf_1.function.subs(x, 0) + ppf_2.function), function_list.append(
            ppf_1.x_max + ppf_2.x_max)
        boundary_list.append(ppf_1.x_max + ppf_2.x_max)
        return UserInterface(output_list, function_list, boundary_list)

    x1, x2 = sp.symbols("x1"), sp.symbols("x2")
    oc_equality = sp.Eq(ppf_1.oc.subs(x, x1), ppf_2.oc.subs(x, x2))
    solutions_of_oc_equality = sp.solve(oc_equality, x1)
    real_oc_solutions = []

    for solution in solutions_of_oc_equality:
        if str(solution).find("I") == -1:
            real_oc_solutions.append(solution)

    x1_solution = real_oc_solutions[-1]

    equation = sp.Eq(x1_solution + x2, x)
    x2_optimum = sp.solve(equation, x2)[0]
    x1_optimum = x1_solution.subs(x2, x2_optimum)

    x1_first = False

    if ppf_1.oc.subs(x, 0) == ppf_2.oc.subs(x, 0):
        equation_1 = sp.Eq(ppf_1.x_max, x1_optimum)
        equation_1_solution = sp.solve(equation_1, x)[-1]

        equation_2 = sp.Eq(ppf_2.x_max, x2_optimum)
        equation_2_solution = sp.solve(equation_2, x)[-1]

        if equation_1_solution < equation_2_solution:
            boundary_point = equation_1_solution
            translation = x1_optimum.subs(x, equation_1_solution)
            piece_2 = ppf_2.function.subs(x, x - translation)
        else:
            boundary_point = equation_2_solution
            translation = x2_optimum.subs(x, equation_2_solution)
            piece_2 = ppf_1.function.subs(x, x - translation)

        boundary_point = boundary_point if type(boundary_point) == sp.Mul else round(float(boundary_point), 3)
        piece_1 = ppf_1.function.subs(x, x1_optimum) + ppf_2.function.subs(x, x2_optimum)
        boundary_point2 = ppf_1.x_max + ppf_2.x_max
        boundary_point2 = boundary_point2 if type(boundary_point2) == sp.Mul else round(float(boundary_point2), 3)

        output_list.append(f"y = {piece_1}, 0 <= x <= {boundary_point}")
        output_list.append(f"y = {piece_2}, {boundary_point} <= x <= {boundary_point2}")
        function_list.append(piece_1), function_list.append(piece_2)
        boundary_list.append(boundary_point), boundary_list.append(boundary_point2)

        return UserInterface(output_list, function_list, boundary_list)

    elif ppf_1.oc.subs(x, 0) < ppf_2.oc.subs(x, 0):
        first_ppf = ppf_1.function
        second_ppf = ppf_2.function
        first_oc = ppf_1.oc
        second_oc = ppf_2.oc
        x1_first = True

    else:
        first_ppf = ppf_2.function
        second_ppf = ppf_1.function
        first_oc = ppf_2.oc
        second_oc = ppf_1.oc
    equation = sp.Eq(first_oc, second_oc.subs(x, 0))
    boundary_point_1 = [root for root in sp.solve(equation, x) if str(root).find("I") == -1 and root >= 0][0]

    # x1_bound and x2_bound are how much of x1 and x2 are exhausted

    if x1_first:
        x1_bound = boundary_point_1
        x2_bound = 0

    else:
        x1_bound = 0
        x2_bound = boundary_point_1

    equation_1 = sp.Eq(ppf_1.x_max - x1_bound, x1_optimum)
    equation_1_solution = sp.solve(equation_1, x)[-1]

    equation_2 = sp.Eq(ppf_2.x_max - x2_bound, x2_optimum)
    equation_2_solution = sp.solve(equation_2, x)[-1]

    boundary_point_2 = equation_1_solution + x1_bound + x2_bound

    piece_2 = ppf_1.function.subs(x, x1_optimum) + ppf_2.function.subs(x, x2_optimum)
    piece_2 = piece_2.subs(x, x - boundary_point_1)

    if equation_1_solution < equation_2_solution:
        piece_3 = ppf_2.function.subs(x, x - boundary_point_2)

    else:
        piece_3 = ppf_1.function.subs(x, x - boundary_point_2)

    piece_1 = first_ppf + second_ppf.subs(x, 0)
    x_max = ppf_1.x_max + ppf_2.x_max
    boundary_point_1 = boundary_point_1 if type(boundary_point_1) == sp.Mul else round(float(boundary_point_1), 3)
    boundary_point_2 = boundary_point_1 if type(boundary_point_2) == sp.Mul else round(float(boundary_point_2), 3)
    x_max = x_max if type(x_max) == sp.Mul else round(float(x_max), 3)

    output_list.append(f"y = {piece_1}, 0 <= x <= {boundary_point_1}")
    output_list.append(f"y = {piece_2}, {boundary_point_1} <= x <= {boundary_point_2}")
    output_list.append((f"y = {piece_3}, {boundary_point_2} <= x <= {x_max}"))

    function_list.append(piece_1), function_list.append(piece_2), function_list.append(piece_3)
    boundary_list.append(boundary_point_1), boundary_list.append(boundary_point_2), boundary_list.append(x_max)


# Non-linear + linear, non-linear is convex
def LinearNonLinear_2(non_linear, linear, x):
    a = sp.symbols("a")

    linear_of_a = linear.function.subs(x, a)
    non_linear_translated = non_linear.function.subs(x, x - a) + linear_of_a

    eq_1 = sp.Eq(x - a, non_linear.x_max)
    non_linear_a_boundary = sp.solve(eq_1, a)[-1]

    eq_2 = sp.Eq(0, non_linear_a_boundary)
    eq_3 = sp.Eq(x, linear.x_max)
    max_changing_point = float(sp.solve(eq_2, x)[-1])
    min_changing_point = float(sp.solve(eq_3, x)[-1])

    if min_changing_point < max_changing_point:
        bound_1 = min_changing_point
        bound_2 = max_changing_point
    else:
        bound_2 = min_changing_point
        bound_1 = max_changing_point

    def get_interval_values(start_point, end_point, mid_point):
        y = sp.sympify("x")
        # since non_linear_x_boundary is a decreasing function, we could check only the start point
        # unless they are equal, then we also need to check the midpoint boundary
        if non_linear_a_boundary.subs(x, start_point) < 0:
            max = 0
        elif non_linear_a_boundary.subs(x, start_point) > 0:
            max = non_linear_a_boundary
        else:
            if non_linear_a_boundary.subs(x, mid_point) > 0:
                max = non_linear_a_boundary
            else:
                max = 0

        # the same for the minimum, but with incr function
        if y.subs(x, end_point) < linear.x_max:
            min = y
        elif y.subs(x, end_point) < linear.x_max:
            min = linear.x_max
        else:
            if y.subs(x, mid_point) < linear.x_max:
                min = y
            else:
                min = linear.x_max

        return max, min

    max_1, min_1 = get_interval_values(0, bound_1, bound_1 / 2)
    max_2, min_2 = get_interval_values(bound_1, bound_2, (bound_1 + bound_2) / 2)
    max_3, min_3 = get_interval_values(bound_2, bound_2 + bound_1, (bound_2 + bound_1 + bound_2) / 2)
    V = linear.function.subs(x, a) + non_linear.function.subs(x, x - a)
    V1max, V1min = V.subs(a, max_1), V.subs(a, min_1)
    V2max, V2min = V.subs(a, max_2), V.subs(a, min_2)
    V3max, V3min = V.subs(a, max_3), V.subs(a, min_3)

    def GetFunction(Vmax, Vmin, bound1, bound2, boundary_list, function_list):
        output_list = []
        equation = sp.Eq(Vmax, Vmin)
        solutions = sp.solve(equation, x)
        solutions = [(solution if type(solution) == sp.Mul else round(float(solution), 3)) for solution in solutions]
        # since one graph could cross another 0, 1 or 2 times only:
        if not solutions:
            midpoint = (bound1 + bound2) / 2
            if Vmax.subs(x, midpoint) >= Vmin.subs(x, midpoint):
                output_list.append(f"y = {Vmax}, {bound1} <= x <= {bound2}")
                function_list.append(Vmax)
            else:
                output_list.append(f"y = {Vmin}, {bound1} <= x <= {bound2}")
                function_list.append(Vmin)
            boundary_list.append(bound1)
            boundary_list.append(bound2)

        elif len(solutions) == 1:
            midpoint = (bound2 + solutions[0]) / 2 if bound2 != solutions[0] else (bound1 + bound2) / 2
            if Vmax.subs(x, midpoint) >= Vmin.subs(x, midpoint):
                if bound2 != solutions[0]:
                    if bound1 != solutions[0]:
                        output_list.append(f"y = {Vmin}, {bound1} <= x <= {solutions[0]}")
                        function_list.append(Vmin)
                        boundary_list.append(bound1)
                        boundary_list.append(solutions[0])
                    output_list.append(f"y = {Vmax}, {solutions[0]} <= x <= {bound2}")
                    function_list.append(Vmax)
                    boundary_list.append(solutions[0])
                    boundary_list.append(bound2)
                else:
                    if bound1 != solutions[0]:
                        output_list.append(f"y = {Vmax}, {bound1} <= x <= {solutions[0]}")
                        function_list.append(Vmax)
                        boundary_list.append(bound1)
                        boundary_list.append(solutions[0])
            else:
                if bound2 != solutions[0]:
                    if bound1 != solutions[0]:
                        output_list.append(f"y = {Vmax}, {bound1} <= x <= {solutions[0]}")
                        function_list.append(Vmax)
                        boundary_list.append(bound1)
                        boundary_list.append(solutions[0])
                    output_list.append(f"y = {Vmin}, {solutions[0]} <= x <= {bound2}")
                    function_list.append(Vmin)
                    boundary_list.append(solutions[0])
                    boundary_list.append(bound2)
                else:
                    if bound1 != solutions[0]:
                        output_list.append(f"y = {Vmin}, {bound1} <= x <= {solutions[0]}")
                        function_list.append(Vmin)
                        boundary_list.append(bound1)
                        boundary_list.append(solutions[0])

        elif len(solutions) == 2:
            midpoint = (solutions[1] + solutions[0]) / 2
            if Vmax.subs(x, midpoint) >= Vmin.subs(x, midpoint):
                if bound1 != solutions[0]:
                    output_list.append(f"y = {Vmin}, {bound1} <= x <= {solutions[0]}")
                    function_list.append(Vmin)
                    boundary_list.append(bound1)
                    boundary_list.append(solutions[0])
                output_list.append(f"y = {Vmax}, {solutions[0]} <= x <= {solutions[1]}")
                function_list.append(Vmax)
                boundary_list.append(solutions[0])
                boundary_list.append(solutions[1])
                if bound1 != solutions[0]:
                    output_list.append(f"y = {Vmin}, {solutions[1]} <= x <= {bound2}")
                    function_list.append(Vmin)
                    boundary_list.append(solutions[1])
                    boundary_list.append(bound2)
            else:
                if bound1 != solutions[0]:
                    output_list.append(f"y = {Vmax}, {bound1} <= x <= {solutions[0]}")
                    function_list.append(Vmax)
                    boundary_list.append(bound1)
                    boundary_list.append(solutions[0])
                output_list.append(f"y = {Vmin}, {solutions[0]} <= x <= {solutions[1]}")
                function_list.append(Vmin)
                boundary_list.append(solutions[0])
                boundary_list.append(solutions[1])
                if bound1 != solutions[1]:
                    output_list.append(f"y = {Vmax}, {solutions[1]} <= x <= {bound2}")
                    function_list.append(Vmax)
                    boundary_list.append(solutions[1])
                    boundary_list.append(bound2)

        return output_list, boundary_list, function_list

    boundary_list = []
    function_list = []
    output_list_1, boundary_list, function_list = GetFunction(V1max, V1min, 0, bound_1, boundary_list, function_list)
    output_list_2, boundary_list, function_list = GetFunction(V2max, V2min, bound_1, bound_2, boundary_list,
                                                              function_list)
    output_list_3, boundary_list, function_list = GetFunction(V3max, V3min, bound_2, bound_1 + bound_2, boundary_list,
                                                              function_list)

    # to merge similar functions and boundaries
    for k in range(len(function_list) - 2, -1, -1):
        if sp.sympify(function_list[k] - function_list[k + 1]) == 0:
            function_list.pop(k + 1)
            boundary_list.pop(2 * k + 1)
            boundary_list.pop(2 * k + 1)

    print(f"""
Let y = {linear.function} be a base
Then, a general vector defining y = {linear.function} is (a, {linear_of_a}),
                                for a within {0, linear.x_max}

 By translate the non-linear function by a general vector we get:

V(x) = max({linear_of_a} + {non_linear.function.subs(x, x - a)})
                       or
V(x) = max({non_linear_translated})

Let's now think of X as an overall amount of X which we allocated between two fields
So, to a linear one, we allocate X = a. Respectively, we allocate X = x - a to the non-linear one.
To satisfy this allocation:

0 <= a <= {linear.x_max}
0 <= x - a <= {non_linear.x_max}, leading to
a >= {non_linear_a_boundary}

Let A(x) be a family of all of the non-linear convex PPF translations. Then,
A(x) = (a: max(0, {non_linear_a_boundary}) <= a <= min(x, {linear.x_max}))

As such, V(x) = max({non_linear_translated}) within A(x) --> max

Because the maximisation is over a closed interval, according to the External Value Theorem,
the maximum or the minimum point will be either at the endpoints or at the interior point (when dy/dx = 0).
However, since V(x) = max(Linear(a) + Non-Linear(x-a)), where Non-Linear is convex, even if an interior point
exists, it will yield a minimum solution, but not a maximum. As such, our optimum solution would lie in midpoints.

Let's look at which x max(0, {non_linear_a_boundary}) and min(x, {linear.x_max}) changes.

0 = {non_linear_a_boundary} leads to x = {max_changing_point}, while
x = {linear.x_max} leads to x = {min_changing_point}

As such, we result in three ranges:
0 <= x <= {bound_1}
{bound_1} < x <= {bound_2}
{bound_2} < x <= {bound_1 + bound_2}

Range 1: 0 <= x <= {bound_1}:
For this region, max(0, {non_linear_a_boundary}) = {max_1}
and min(x, {linear.x_max}) = {min_1}

V(x) = {V}

a = {max_1}:
V1(max) = {V.subs(a, max_1)}

a = {min_1}:
V1(min) = {V.subs(a, min_1)}

V1(max) > V1(min):

{V1max} > {V1min}

On a range 0 <= x <= {bound_1} gives us:
\n{'\n'.join(str(item) for item in output_list_1)}

Range 2: {bound_1} <= x <= {bound_2}:
For this region, max(0, {non_linear_a_boundary}) = {max_2}
and min(x, {linear.x_max}) = {min_2}

V(x) = {V}

a = {max_2}:
V2(max) = {V.subs(a, max_2)}

a = {min_2}:
V2(min) = {V.subs(a, min_2)}

V2(max) > V2(min):

{V2max} > {V2min}

On a range {bound_1} <= x <= {bound_2} gives us:
\n{'\n'.join(str(item) for item in output_list_2)}

Range 3: {bound_2} <= x <= {bound_1 + bound_2}:
For this region, max(0, {non_linear_a_boundary}) = {max_3}
and min(x, {linear.x_max}) = {min_3}
V(x) = {V}

a = {max_3}:
V2(max) = {V.subs(a, max_3)}

a = {min_3}:
V2(min) = {V.subs(a, min_3)}

V3(max) > V3(min):

{V3max} > {V3min}

On a range {bound_2} <= x <= {bound_1 + bound_2} gives us:
\n{'\n'.join(str(item) for item in output_list_3)}

Thus, we are left with the following functions:
\n{'\n'.join(str(item) for item in output_list_1)} \n{'\n'.join(str(item) for item in output_list_2)} \n{'\n'.join(str(item) for item in output_list_3)}

By combining them, we get:
""")
    output_list = []
    for i in range(len(function_list)):
        output_list.append(f"y = {function_list[i]}, {boundary_list[2 * i]} <= x <= {boundary_list[2 * i + 1]}")

    return UserInterface(output_list, function_list, boundary_list)


# Two linears
def Linear(ppf_1, ppf_2, x):
    function_list = []
    output_list = []
    boundary_list = []
    y_max = ppf_1.y_max + ppf_2.y_max
    x_max = ppf_1.x_max + ppf_2.x_max

    if ppf_1.oc > ppf_2.oc:
        piece_1 = sp.sympify(y_max - ppf_2.oc * x)
        boundary_point = ppf_2.x_max
        piece_2 = ppf_1.function.subs(x, x - boundary_point)

        output_list.append(f"y = {piece_1}, 0 <= x <= {boundary_point}")
        output_list.append(f"y = {piece_2}, {boundary_point} <= x <= {x_max}")
        function_list.append(piece_1), function_list.append(piece_2)
        boundary_list.append(boundary_point), boundary_list.append(x_max)

    elif ppf_1.oc < ppf_2.oc:
        piece_1 = sp.sympify(y_max - ppf_1.oc * x)
        boundary_point = ppf_1.x_max
        piece_2 = ppf_2.function.subs(x, x - boundary_point)

        output_list.append(f"y = {piece_1}, 0 <= x <= {boundary_point}")
        output_list.append(f"y = {piece_2}, {boundary_point} <= x <= {x_max}")
        function_list.append(piece_1), function_list.append(piece_2)
        boundary_list.append(boundary_point), boundary_list.append(x_max)

    else:
        final_function = sp.sympify(y_max - ppf_1.oc * x)
        output_list.append(f"y = {final_function}, 0 <= x <= {x_max}")
        function_list.append(final_function)
        boundary_list.append(x_max)

    return UserInterface(output_list, function_list, boundary_list)


# Two non-linears, both convex
def NonLinear_3(ppf_1, ppf_2, x):
    print("The combination is too complex to compute for now")


# Two non-linears, one convex
def NonLinear_2(ppf_1, ppf_2, x):
    print("The combination is too complex to compute for now")


Main()
