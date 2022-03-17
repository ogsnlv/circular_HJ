
namespace circular_HJ
{
    class Program
    {
        private const double Bo = 1;
        private const int N = 100000, n = 4;
        static readonly double step = (2 - 50) / Convert.ToDouble(N);

        static double DerH(double r)
        {
            return (H(r + step) - H(r)) / step;
        }

        static double DerU(double r)
        {
            return (U(r + step) - U(r)) / step;
        }

        static double Derh1(double Fr, double r)
        {
            return (h1(Fr, r + step) - h1(Fr, r)) / step;
        }

        static double Fh(double h, double u, double w, double I, double r)
        {
            return 1 / (1 - U(r) * U(r) / H(r)) * ((I + DerU(r)) * (U(r) / H(r) * h - u) + U(r) / H(r) * (1 / r * (u * H(r) + h * U(r)) + u * DerH(r) + H(r) * n * w / r));
        }

        static double Fr(double r)
        {
            return Math.Pow(U(r), 2) / H(r);
        }

        static double Fu(double h, double u, double w, double I, double r)
        {
            return 1 / (U(r) * U(r) - H(r)) * ((I + DerU(r)) * (h - u * U(r)) + 1 / r * (u * H(r) + h * U(r)) + u * DerH(r) + H(r) * n * w / r);
        }

        static double Fw(double h, double u, double w, double I, double r)
        {
            return 1 / U(r) * (n * h - w * (u / r + I));
        }

        static double H(double r)
        {
            return 1 - 1 / (Math.Pow(r, 2) - 2);
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

            double alpha, det;
            double[] v1 = { 1, 0, 0 }, v2 = { 0, 0, 1 }, v3 = new double[3];

            for (int i = N - 1; i > 0; i--)
            {

                alpha = (1 + 2 * Bo * U(r[i]) + Math.Sqrt(Math.Pow(1 + 2 * Bo * U(r[i]), 2)) + 8 * Fr(r[i])) / (4 * Fr(r[i]));

                RungeKutt(ref v1, I, r[i]);
                RungeKutt(ref v2, I, r[i]);

                Normalization(ref v1);
                Normalization(ref v2);

                Orthogonalization(v1, ref v2);

                Normalization(ref v1);
                Normalization(ref v2);

// jump conditions
                if (r[i] < 10)
                {
                    v3[0] = 1 / (H(r[i]) - Math.Pow(U(r[i]), 2) + Bo / r[i]) * (-Bo / r[i] * (DerH(r[i]) -
                        Derh1(Fr(r[i]), r[i])) - H(r[i]) * (1 - 1 / alpha) * (Bo / Math.Pow(r[i], 2) * (Math.Pow(n, 2) - 1) - 2 * U(r[i]) * I));

                    v3[1] = 1 / (U(r[i]) - H(r[i]) / U(r[i]) - Bo / (r[i] * U(r[i]))) * (-Bo / (r[i] * H(r[i])) * (DerH(r[i]) - Derh1(Fr(r[i]), r[i])) - (1 - 1 / alpha) * (Bo / (r[i] * U(r[i])) +
                        Bo / Math.Pow(r[i], 2) * (Math.Pow(n, 2) - 1) + U(r[i]) * I + I * H(r[i]) / U(r[i])));
                    v3[2] = -n * U(r[i]) * (alpha - 1);

                    Normalization(ref v3);

                    det = v1[0] * v2[1] * v3[2] + v1[1] * v2[2] * v3[0] + v1[2] * v2[0] * v3[1] - v1[2] * v2[1] * v3[0] - v1[1] * v2[0] * v3[2] - v1[0] * v2[2] * v3[1];

                    Console.WriteLine("{0};{1}", r[i], det);
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

            k1 = Fh(h, u, w, I, r);
            k2 = Fh(h + step * 0.5 * k1, u + step * 0.5 * k1, w + step * 0.5 * k1, I, r + step * 0.5 * k1);
            k3 = Fh(h + (step * 0.5 * k2), u + (step * 0.5 * k2), w + (step * 0.5 * k2), I, r + step * 0.5 * k2);
            k4 = Fh(h + step * k3, u + step * k3, w + step * k3, I, r + step * k3);

            l1 = Fu(h, u, w, I, r);
            l2 = Fu(h + step * 0.5 * l1, u + step * 0.5 * l1, w + step * 0.5 * l1, I, r + step * 0.5 * k1);
            l3 = Fu(h + step * 0.5 * l2, u + step * 0.5 * l2, w + step * 0.5 * l2, I, r + step * 0.5 * k2);
            l4 = Fu(h + step * l3, u + step * l3, w + step * l3, I, r + step * k3);

            m1 = Fw(h, u, w, I, r);
            m2 = Fw(h + step * 0.5 * m1, u + step * 0.5 * m1, w + step * 0.5 * m1, I, r + step * 0.5 * k1);
            m3 = Fw(h + step * 0.5 * m2, u + step * 0.5 * m2, w + step * 0.5 * m2, I, r + step * 0.5 * k2);
            m4 = Fw(h + step * m3, u + step * m3, w + step * m3, I, r + step * k3);

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

        static double h1(double Fr, double r)
        {
            return 4 * Fr * H(r) / (1 + 2 * Bo * U(r) + Math.Sqrt(Math.Pow(1 + 2 * Bo * U(r), 2) + 8 * Fr));
        }
    }
}