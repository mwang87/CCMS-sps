// Module Includes
#include "ExecBase.h"

// System Includes
#include <stdio.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
ExecBase::ExecBase(void) :
  m_isValid(false), m_name("ExecBase")
{
}

// -------------------------------------------------------------------------
ExecBase::ExecBase(const ParameterList & params) :
  m_params(params), m_name("ExecBase"), m_isValid(false)
{
}

// -------------------------------------------------------------------------
ExecBase::~ExecBase(void)
{
}

// -------------------------------------------------------------------------
bool ExecBase::isValid(void) const
{
  return m_isValid;
}

// -------------------------------------------------------------------------
string ExecBase::getName(void) const
{
  return m_name;
}

// -------------------------------------------------------------------------
string ExecBase::getType(void) const
{
  return m_type;
}

// -------------------------------------------------------------------------
void ExecBase::setName(std::string name)
{
  m_name = name;
}

// -------------------------------------------------------------------------
std::string ExecBase::makeName(std::string baseFilename,
                               int numSplit) const
{
  std::string filename(baseFilename);
  char buf[256];
  sprintf(buf, "%d", numSplit + 1);

  filename += "_";
  filename += buf;

  return filename;
}

