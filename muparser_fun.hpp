#include <muParser.h>

#include <memory>
#include <string>

class MuparserFun
{
public:
  MuparserFun(const MuparserFun &m)
    : m_parser(m.m_parser)
  {
    m_parser.DefineVar("x", &m_var);
    m_parser.DefineVar("y",&m_var2);
  };

  MuparserFun(const std::string &s)
  {
    try
      {
        m_parser.DefineVar("x", &m_var);
        m_parser.DefineVar("y", &m_var2);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double operator()(const double &x,const double &y)
  {
    m_var = x;
    m_var2=y;
    double z = m_parser.Eval();
    return z;
  };

private:
  double     m_var;
  double m_var2;
  mu::Parser m_parser;
};