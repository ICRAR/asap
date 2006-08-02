#ifndef ASAPLOGGER_H
#define ASAPLOGGER_H

#include <string>
#include <iostream>

namespace asap {
/**
  * This class provides the logging within asap. All other classes which need to log
  * should inherit from this.
  * @brief The ASAP logging class
  * @author Malte Marquarding
  * @date $Date$
  * @version
  */
class Logger {
public:
  /**
   * Default Constructor
   **/
  Logger();
  
  /**
   * Constructor with switch to enable/disable logging
   * @param[in] enabled indicating the deafult state
   */  
  Logger(bool enabled);

  /*
   * Destructor
   */
  virtual ~Logger();
  /**
   * push another message into the logger
   * @param[in] the message
   * @param[in] whether to add a newline character at the end
   */  
  void pushLog(const std::string& s, bool newline=true) const;
  /**
   * pop the message form the logger
   * @returns the log message string
   */  
  std::string popLog() const;
  /**
   * enable logging
   */
  virtual void enableLog();
  /**
   * disable logging
   */
  virtual void disableLog();
  
private:
  static std::string log_;
  bool enabled_;
};

} // namespace

#endif

