clc;
clear;

// =============================
// 1. Gaussian Model Function
// =============================
function f = gaussian_model(x, a)
    A     = a(1);
    mu    = a(2);
    sigma = max(a(3), 0.1);   // prevent sigma collapse

    f = A * exp(-((x - mu).^2) ./ (2*sigma^2));
endfunction


// =============================
// 2. Jacobian Function (J = dr/da)
// =============================
function J = gaussian_jacobian(x, a)
    A     = a(1);
    mu    = a(2);
    sigma = max(a(3), 0.1);   

    exp_term = exp(-((x - mu).^2) ./ (2*sigma^2));

    // ∂f/∂A
    J1 = exp_term;

    // ∂f/∂mu
    J2 = A * exp_term .* ((x - mu) / (sigma^2));

    // ∂f/∂sigma
    J3 = A * exp_term .* ((x - mu).^2 / (sigma^3));

    //Jacobian of residual r = y - f → negative
    J = -[J1, J2, J3];
endfunction


// =============================
// 3. Levenberg–Marquardt Solver
// =============================
function a = gauss_newton(x, y, a0, max_iter, tol)

    a = a0;
    lambda = 1;   // better start

    for k = 1:max_iter

        // Current residual
        f = gaussian_model(x, a);
        r = y - f;
        cost = r' * r;

        // Jacobian
        J = gaussian_jacobian(x, a);
        H = J' * J;
        g = J' * r;

        // LM step
        delta = (H + lambda*eye(3)) \ g;

        // Trial update (NO alpha)
        a_trial = a + delta;

        // Clamp
        a_trial(3) = max(a_trial(3), 0.1);
        a_trial(2) = max(min(a_trial(2), max(x)), min(x));

        // New cost
        r_new = y - gaussian_model(x, a_trial);
        cost_new = r_new' * r_new;

        // Decision rule
        if cost_new < cost then
            a = a_trial;
            lambda = lambda / 10;   // trust Gauss-Newton more
        else
            lambda = lambda * 10;   // trust gradient descent more
        end

        // Stop condition
        if norm(delta) < tol then
            break;
        end

    end

endfunction


// =============================
// 4. Main Script
// =============================

// Data

x = [-5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]';

y = [0.81, -1.24, 2.15, 0.33, -0.55, 1.9, -0.12, 1.44, 0.98, -1.1, 2.5, 0.67, 1.12, -0.88, 0.45, 1.89, 2.33, 3.41, 4.12, 4.88, 6.12, 7.9, 9.44, 11.55, 14.22, 17.88, 21.33, 25.12, 30.88, 34.41, 39.55, 45.12, 52.33, 58.9, 63.44, 70.12, 75.88, 82.33, 87.41, 91.55, 94.12, 96.88, 98.33, 99.12, 101.44, 102.12, 101.88, 99.55, 100.9, 101.33, 102.45, 100.22, 97.45, 94.12, 91.88, 86.41, 81.33, 76.54, 70.12, 64.88, 59.33, 53.41, 47.9, 41.22, 35.88, 31.44, 26.5, 22.12, 17.88, 15.33, 11.41, 9.55, 7.33, 6.12, 5.88, 4.12, 2.9, 1.55, 2.12, 0.88, -1.22, 1.45, 0.12, -0.55, 1.88, 0.44, -1.12, 0.9, 1.33, -0.12, 1.55, -0.88, 0.45, 1.22, -1.44, 0.33, 1.1, -0.21, 0.88, -1.33, 0.55]';

// Initial guess
a0 = [105; -0.4; 1.2];

// Settings
max_iter = 100;
tol = 1e-6;

// Solve
a = gauss_newton(x, y, a0, max_iter, tol);

// Display result
disp("Final parameters [A, mu, sigma]:");
disp(a);


// =============================
// 5. Plot
// =============================
x_fit = linspace(min(x), max(x), 200)';
y_fit = gaussian_model(x_fit, a);

plot(x, y, 'bo', x_fit, y_fit, '-r');
xlabel("x");
ylabel("y");
title("Gaussian Fit )");
legend("Data points", "Fitted Gaussian");
