using MathNet.Numerics.LinearAlgebra.Double;

namespace circular_HJ
{
    class Program
    {
        private const double Bo = 1;
        private const int N = 100000, n = 4;
        static readonly double step = (2 - 50) / Convert.ToDouble(N);
        static double[,] a = new double[3,3];
        static double[] b = new double[3];
        static double c_1_left = 1.26, c_2_left = 2.28, c_1_right = 1.2, c_2_right = 3;

        static double Froude(double r)
        {
            return Math.Pow(U(r), 2) / H(r);
        }

        static double H(double r)
        {
            return 1 - 1 / (Math.Pow(r, 2) - 2);
        }

        static double H_derivative(double r)
        {
            return (H(r + step) - H(r)) / step;
        }

        static double Kramer(double[,] a, double[] b, int counter)
        {
            Matrix m = DenseMatrix.OfArray(a);
            double[] x = new double[3];
            double det = m.Determinant();
            for (int i = 0; i < 3; i++)
            {
                a[i, counter] = b[counter];
            }
            m = DenseMatrix.OfArray(a);
            double det_i = m.Determinant();
            return det_i / det;

        }
        static void Main(string[] args)
        {
            int I = 0;
            double[] r = new double[N];
            r[N - 1] = 50;
            for (int i = N - 1; i > 0; i--)
            {
                r[i - 1] = r[i] + step;
            }

            double alpha;
            double[] v1 = { 1, 0, 0 }, v2 = { 0, 0, 1 }, v3 = new double[3];

            for (int i = N - 1; i > 0; i--)
            {

                alpha = (1 + 2 * Bo * U(r[i]) + Math.Sqrt(Math.Pow(1 + 2 * Bo * U(r[i]), 2)) + 8 * Froude(r[i])) / (4 * Froude(r[i]));

                RungeKutt(ref v1, I, r[i]);
                RungeKutt(ref v2, I, r[i]);

                Normalization(ref v1);
                Normalization(ref v2);

                Orthogonalization(v1, ref v2);

                Normalization(ref v1);
                Normalization(ref v2);

                // jump conditions
                Matrix m;
                if (r[i] < 10)
                {
                    // v3[0] = 1 / (H(r[i]) - Math.Pow(U(r[i]), 2) + Bo / r[i]) * (-Bo / r[i] * (H_derivative(r[i]) - h1_derivative(Froude(r[i]), r[i])) - H(r[i]) * (1 - 1 / alpha) * (Bo / Math.Pow(r[i], 2) * (Math.Pow(n, 2) - 1) - 2 * U(r[i]) * I));

                    // v3[1] = 1 / (U(r[i]) - H(r[i]) / U(r[i]) - Bo / (r[i] * U(r[i]))) * (-Bo / (r[i] * H(r[i])) * (H_derivative(r[i]) - h1_derivative(Froude(r[i]), r[i])) - (1 - 1 / alpha) * (Bo / (r[i] * U(r[i])) + Bo / Math.Pow(r[i], 2) * (Math.Pow(n, 2) - 1) + U(r[i]) * I + I * H(r[i]) / U(r[i])));
                    //v3[2] = -n * U(r[i]) * (alpha - 1);

                    Normalization(ref v3);
                    // det = v1[0] * v2[1] * v3[2] + v1[1] * v2[2] * v3[0] + v1[2] * v2[0] * v3[1] - v1[2] * v2[1] * v3[0] - v1[1] * v2[0] * v3[2] - v1[0] * v2[2] * v3[1];
                    m = DenseMatrix.OfRowArrays(v1, v2, v3);
                    Console.WriteLine("{0};{1}", r[i], m.Determinant());
                }
            }
        }

        static void Normalization(ref double[] v)
        {
            double norm = Math.Sqrt(Math.Pow(v[0], 2) + Math.Pow(v[1], 2) + Math.Pow(v[2], 2));

            for (int i = 0; i < 3; i++)
            {
                v[i] = v[i] / norm;
            }

        }

        static void Orthogonalization(double[] v1, ref double[] v2)
        {
            double product1 = 0, product2 = 0;

            for (int i = 0; i < 3; i++)
            {
                product1 += v1[i] * v2[i];
                product2 += v1[i] * v1[i];
            }
            for (int i = 0; i < 3; i++)
            {
                v2[i] = v2[i] - product1 / product2 * v1[i];
            }
        }

        static void RungeKutt(ref double[] v, double I, double r)
        {
            double h, u, w;
            double k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4;

            h = v[0];
            u = v[1];
            w = v[2];

            k1 = h_function(h, u, w, I, r);
            k2 = h_function(h + step * 0.5 * k1, u + step * 0.5 * k1, w + step * 0.5 * k1, I, r + step * 0.5 * k1);
            k3 = h_function(h + (step * 0.5 * k2), u + (step * 0.5 * k2), w + (step * 0.5 * k2), I, r + step * 0.5 * k2);
            k4 = h_function(h + step * k3, u + step * k3, w + step * k3, I, r + step * k3);

            l1 = u_function(h, u, w, I, r);
            l2 = u_function(h + step * 0.5 * l1, u + step * 0.5 * l1, w + step * 0.5 * l1, I, r + step * 0.5 * k1);
            l3 = u_function(h + step * 0.5 * l2, u + step * 0.5 * l2, w + step * 0.5 * l2, I, r + step * 0.5 * k2);
            l4 = u_function(h + step * l3, u + step * l3, w + step * l3, I, r + step * k3);

            m1 = w_function(h, u, w, I, r);
            m2 = w_function(h + step * 0.5 * m1, u + step * 0.5 * m1, w + step * 0.5 * m1, I, r + step * 0.5 * k1);
            m3 = w_function(h + step * 0.5 * m2, u + step * 0.5 * m2, w + step * 0.5 * m2, I, r + step * 0.5 * k2);
            m4 = w_function(h + step * m3, u + step * m3, w + step * m3, I, r + step * k3);

            h = h + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            u = u + step * (l1 + 2 * l2 + 2 * l3 + l4) / 6;
            w = w + step * (m1 + 2 * m2 + 2 * m3 + m4) / 6;

            v[0] = h;
            v[1] = u;
            v[2] = w;

            Normalization(ref v);

        }
        static double U(double r)
        {
            return 1 / Math.Pow(r, 1.1) * (1 + 1 / (Math.Pow(r, 2) - 2));
        }

        static double U_derivative(double r)
        {
            return (U(r + step) - U(r)) / step;
        }

        static double h1(double Froude, double r)
        {
            return 4 * Froude * H(r) / (1 + 2 * Bo * U(r) + Math.Sqrt(Math.Pow(1 + 2 * Bo * U(r), 2) + 8 * Froude));
        }

        static double h1_derivative(double Fr, double r)
        {
            return (h1(Fr, r + step) - h1(Fr, r)) / step;
        }
        static double h_function(double h, double u, double w, double I, double r)
        {
            return Kramer(a, b, 0);
            // return 1 / (1 - U(r) * U(r) / H(r)) * ((I + U_derivative(r)) * (U(r) / H(r) * h - u) + U(r) / H(r) * (1 / r * (u * H(r) + h * U(r)) + u * H_derivative(r) + H(r) * n * w / r));
        }

        static double u_function(double h, double u, double w, double I, double r)
        {
            return Kramer(a, b, 1);
            //return 1 / (U(r) * U(r) - H(r)) * ((I + U_derivative(r)) * (h - u * U(r)) + 1 / r * (u * H(r) + h * U(r)) + u * H_derivative(r) + H(r) * n * w / r);
        }

        static double w_function(double h, double u, double w, double I, double r)
        {
            return Kramer(a, b, 2);
            //return 1 / U(r) * (n * h - w * (u / r + I));
        }
    }

    class StaticJump
    {
        const int N = (int)1e4;

        static double c_1_left = 1.26, c_2_left = 2.28, c_1_right = 1.2, c_2_right = 3;

        static double Forward(double r, double h)
        {
            return (h - (c_2_left * Math.Pow(r, 2))) / (r * ((Math.Pow(h, 3) * Math.Pow(r, 2)) - c_1_left));
        }

        static void Main()
        {

            int choice = 0, counter = 0;

            double
                difference_1, difference_2,
                We, Re,
                a_s = 0.002,
                a = Math.Pow(a_s, 1.25) * Math.Pow(0.025, -0.25),
                nu = 32e-6,
                g = 9.81,
                rho = 1120, Bo = 0.76,
                sigma = rho * g * Math.Pow(a_s, 2) / Bo,
                temp = 0,
                B_j;

            double Q, H_c, R_c, U_c, Bo_1, H, r_0, step;

            We = Math.Pow(10, 1.13101);
            Re = Math.Pow(10, 2.60361);
            double[] r = new double[N], h_left = new double[N], u_left = new double[N], h_right = new double[N], u_right = new double[N];

            Q = Re * nu * a;
            H_c = Math.Pow(Q * nu / (2 * Math.PI * g), 0.25);
            U_c = Math.Pow(Q * nu * Math.Pow(g, 2) / (2 * Math.PI), 0.125);
            R_c = Math.Pow(Q / (2 * Math.PI), 0.625) / Math.Pow(g * Math.Pow(nu, 3), 0.125);

            Bo_1 = sigma / (Math.Pow(U_c, 2) * R_c * rho);
            //H = Math.Pow(rho * Q * Q / (sigma * We), 0.333);
            H = 0.0122;
            r_0 = 0.3155 * a * Math.Pow(Re, 0.33) / R_c;
            step = (10 - r_0) / (double)N;

            r[0] = r_0;
            h_left[0] = 5.74 * Math.Pow(r[0], 2);
            u_left[0] = 1 / (h_left[0] * r[0]);
            h_right[N - 1] = H / H_c;
            u_right[N - 1] = 1 / (h_right[N - 1] * r[N - 1]);

            for (int i = 0; i < N - 1; i++)
            {
                r[i + 1] = r[i] + step;
                h_left[i + 1] = RG(r[i], h_left[i], step, true);
                u_left[i + 1] = 1 / (h_left[i + 1] * r[i + 1]);
                // double I=0; if (Math.Abs(h_left[i + 1] - h_left[i]) / step > 10) { I = i + 1; break; }
            }

            for (int i = N - 2; i >= 0; i--)
            {
                h_right[i] = RG(r[i + 1], h_right[i + 1], -step, false);
                u_right[i] = 1 / (h_right[i] * r[i]);

            }

            choice = 0;
            switch (choice)
            {
                case 0:
                    {
                        counter = 0;
                        for (int i = 0; i < N; i++)
                        {
                            difference_1 = Math.Pow(u_left[i], 2) * c_1_left * h_left[i] - (Math.Pow(u_right[i], 2) * c_1_right * h_right[i]) + 0.5 * (Math.Pow(h_left[i], 2) - Math.Pow(h_right[i], 2)) - Bo_1 * (h_left[i] - h_right[i]) / r[i];
                            Console.WriteLine("{0};{1}", r[i], difference_1);

                            if (difference_1 * temp < 0) break;

                            temp = difference_1;
                            counter = i;
                        }

                        //Console.WriteLine("no but maybe:{0};{1};{2}", Math.Log10(We), Math.Log10(Re), r[counter]);
                        break;
                    }

                case 1:
                    {
                        for (int i = 0; i < N; i++)
                        {
                            //r[i] = 0.002 + i * 0.001;
                            B_j = rho * g * r[i] * (h_left[i] - h_right[i]) / sigma;

                            //difference_2 = r[i] * g * Math.Pow(H * a / Q, 2) * (1 + 2 / B_j) + Math.Pow(a / Math.PI, 2) / (2 * r[i] * H) - 0.10132 + 0.1297 * Math.Pow(r[i] / a, 1.5) * Math.Pow(Re, -0.5);
                            difference_2 = r[i] * g * Math.Pow(H * a / Q, 2) * (1 + 2 / B_j) + Math.Pow(a / Math.PI, 2) / (2 * r[i] * H) - 0.01676 * Math.Pow(Math.Pow(r[i] / a, 3) * Math.Pow(Re, -1) + 0.1826, -1);
                            Console.WriteLine("{0};{1}", r[i], difference_2);
                            if (difference_2 * temp < 0) break;

                            temp = difference_2;

                        }
                        break;
                    }
            }
        }

        static double RG(double r, double h, double step, bool flag)
        {
            double k_1, k_2, k_3, k_4;
            if (flag == true)
            {
                k_1 = Forward(r, h);
                k_2 = Forward(r + 0.5 * step, h + 0.5 * step * k_1);
                k_3 = Forward(r + 0.5 * step, h + 0.5 * step * k_2);
                k_4 = Forward(r + step, h + step * k_3);
            }

            else
            {
                k_1 = Reverse(r, h);
                k_2 = Reverse(r + 0.5 * step, h + 0.5 * step * k_1);
                k_3 = Reverse(r + 0.5 * step, h + 0.5 * step * k_2);
                k_4 = Reverse(r + step, h + step * k_3);
            }

            return h + 0.166 * step * (k_1 + 2 * (k_2 + k_3) + k_4);
        }

        static double Reverse(double r, double h)
        {
            return (h - c_2_right * Math.Pow(r, 2)) / (r * (Math.Pow(h, 3) * Math.Pow(r, 2) - c_1_right));
        }
    }
}

