#ifndef SDLOG_H
#define SDLOG_H

#include <string>
#include <iostream>

namespace asap {

class SDLog {
public:
  SDLog();
  SDLog(bool enabled);
  void pushLog(const std::string& s, bool newline=true) const;
  //void pushLog(const char* cs);
  std::string popLog() const;
  virtual void enableLog();
  virtual void disableLog();
private:
  static std::string log_;
  bool enabled_;
};

}

#endif

