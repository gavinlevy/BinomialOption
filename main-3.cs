/******************************************************************************

                            Online C# Compiler.
                Code, Compile, Run and Debug C# program online.
Write your code in this editor and press "Run" button to execute it.

*******************************************************************************/

using System;
class Project3 {
  static void Main() {
    
  }
  static void option(){
      
        BinomialOption p = new BinomialOption(100, 1, 0.07,0.02,0.2);
        p.print();
        double s0=88;
        double t0=0;
        double Vtarget=18;
        double Vtol=0.01;
        p.Valuation(s0,t0);
        p.impliedVolatility(s0,t0,Vtarget,Vtol);
    }
}


class BinomialOption
{
  private bool type;
  private bool sort;
  private double K;
  private double T;
  private double rf;
  private double q;
  private double σ;

  public BinomialOption (bool OptionType, bool Optionsort, double Strike, double ExpireTime, double RiskFreeRate, double ContinuousDividendYield, double Volatility)
  {
    type = OptionType;
    sort = Optionsort;
    K = Strike;
    T = ExpireTime;
    rf = RiskFreeRate;
    q = ContinuousDividendYield;
    σ = Volatility;
  }

  public string getOptionType()
  {
    if(type == true) { return "Call Option."; } else { return "Put Option."; }
  }
  public string getOptionsort()
  {
    if(sort == true) { return "American Option."; } else { return "European Option."; }
  }
  public double getStrike ()                  {return K;}
  public double getExpirationTime ()          {return T;}
  public double getRiskFreeRate ()            {return rf;}
  public double getContinuousDividendYield () {return q;}
  public double getVolatility ()              {return σ;}

  public void print ()
  {
    if(type == true)       { Console.WriteLine("This option is a call"); } else if(type == false) { Console.WriteLine("This option is a put"); }
    if(sort == true)       { Console.WriteLine("This option is an American."); } else if(sort == false) { Console.WriteLine("This option is an European."); }
    Console.WriteLine("The Strike Price is " + K);
    Console.WriteLine("The Expiration Time is " + T + " years");
    Console.WriteLine("The Risk Free Rate is " + rf * 100 + "%");
    Console.WriteLine("The Continuous Dividend Yield is " + q * 100 + "%");
    Console.WriteLine("The Volatility is " + σ * 100 + "%");
  }

 // Calculate the parameter values for u, d, the discount factor df and the risk-neutral probabilities p rn, q rn, etc
    public double Valuation(double S0, double t0, int n)
  {
      double dt = (T - t0) / (double)n;
      double u = Math.Exp(σ * Math.Sqrt(dt));
      double d = Math.Exp(- σ * Math.Sqrt(dt));
      double df = Math.Exp(- rf * dt);
      double p_rn = (Math.Exp(dt * (rf - q)) - d) / (u - d);
      double q_rn = (u - Math.Exp(dt * (rf - q))) / (u - d);

 // Memory allocation
      double[][] S = null;
      double[][] V = null;
      S = new double[n + 1][];
      V = new double[n + 1][];
      for (int i = 0; i <= n; i++)
      {
          S[i] = new double[i + 1];
          V[i] = new double[i + 1];
        }

 // Stock prices in nodes
     // Set the values of the stock prices at every node
      S[0][0] = S0;                      // first node
      for (int i = 1; i <= n; i++)       // note loop account
        {
          for (int j = 1; j <= i; j++)   // note loop account
            {
              S[i][j] = S0 * Math.Pow(u, (i + 1 - j)) * Math.Pow(d, j);  // compute the value of S[i][j]
            } 
      }

// Terminal payoff
     //set the terminal payoff at the stepn = n
      for (int j = 0; j <= n; j++)
      {
          if(type == true)       { V[n][j] = Math.Max((S[n][j] - K), 0); }    // terminal payoff depends on value of S[n][j] for Call Option
          else if(type == false) { V[n][j] = Math.Max((K - S[n][j]), 0); }    // terminal payoff depends on value of S[n][j] for Put Option
        }

// Discounted expectation and valuation tests
      for (int i = n - 1; i >= 0; i--)    // work backwards
      {
          for (int j = 0; j <= i; j++)
          {   // discounted expectation
              double Vtmp = df * (p_rn * V[i + 1][j] + q_rn * V[i + 1][j + 1]);
              //overwrite if it's American Option.
              if(sort == true)
              {
                  if(type == true)  { Vtmp = Math.Max(Vtmp, Math.Max((S[i][j] - K), 0)); }
                  if(type == false) { Vtmp = Math.Max(Vtmp, Math.Max((K - S[i][j]), 0)); }
              }
              
              V[i][j] = Vtmp;
          }
        }
 // Return the value of V[0][0] as the fair value of the option
      return V[0][0];
  }
  
// Calculate and return Greeks
  public double Del(double S0, double t0, int n)
  {
      double delta = (Valuation((S0 + S0 / 1000), t0, n) - Valuation((S0 - S0 / 1000), t0, n)) / (2 * S0 / 1000);
      
      return delta;
  }
  
  public double Gam(double S0, double t0, int n)
  {
      double gamma = (Valuation((S0 + S0 / 1000), t0, n) + Valuation((S0 - S0 / 1000), t0, n) - 2 * Valuation(S0, t0, n)) / Math.Pow(S0 / 1000, 2);
      
      return gamma;
  }
  
  public double Veg(double S0, double t0, int n)
  {
      double sigma1 = σ + σ / 1000;
      double sigma2 = σ - σ / 1000;
      
      double σtmp = σ;
      
      σ = sigma1;
      double vvega1 = Valuation(S0, t0, n);
      
      σ = sigma2;
      double vvega2 = Valuation(S0, t0, n);
      
      double vega = (vvega1 - vvega2) / (2 * σtmp / 1000);
      
      σ = σtmp;
      
      return vega;
  }
  
  public double rho(double S0, double t0, int n)
  {
      double r1 = rf + rf / 1000;
      double r2 = rf - rf / 1000;
      
      double rtmp = rf;
      
      rf = r1;
      double vrho1 = Valuation(S0, t0, n);
      
      rf = r2;
      double vrho2 = Valuation(S0, t0, n);
      
      double rho = (vrho1 - vrho2) / (2 * rtmp / 1000);
      
      rf = rtmp;
      
      return rho;
  }
  
  public double Thet(double S0, double t0, int n)
  {
      double theta = (Valuation(S0, (t0 + (T - t0) / (100 * (double)n)), n) - Valuation(S0, t0, n)) / ((T - t0) / (100 * (double)n));
      
      return theta;
  }


// Calculate and return the implied volatility
  public double impliedVolatility(double S0, double t0, double TargetPrice, double tolerance, int n)
  {
      for(int i = 0; i <= n; i++)
      {
          double ivalue = Valuation(S0, t0, n);
          double ivega = Veg(S0, t0, n);
          
          if (Math.Abs(ivalue - TargetPrice) < tolerance) { break; }
          
          σ = σ - (ivalue - TargetPrice) / ivega;
      }
      return (σ * 100);
  }
}
