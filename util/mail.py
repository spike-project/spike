#!/usr/bin/env python 
# encoding: utf-8

from __future__ import print_function
import smtplib
import sys
if sys.version_info[0] < 3:
  from email.MIMEMultipart import MIMEMultipart
  from email.MIMEBase import MIMEBase
  from email.MIMEText import MIMEText
  from email import Encoders
else:
  from email import encoders as Encoders
  from email.mime.base import MIMEBase
  from email.mime.text import MIMEText
  from email.mime.multipart import MIMEMultipart

import os.path as op
import unittest
###

class GMAIL(object):
    '''
    Class for sending mails with Gmail with smtp protocol.
    input:
        gmail_user : name of the Gmail account
        gmail_pwd  : password of the Gmail account
        to: destination of the mail
        subject: subject of the mail
        text: text to be sent
        attach: Attached document
    Usage: 
        gm = GMAIL()
        gm.send(to = 'gmalert67@gmail.com', subject = 'test gmail', text = "hello", attach = None)
    '''
    def __init__(self, gmail_user = "gmalert67@gmail.com", gmail_pwd = "igbmcalert" ):
        self.gmail_user = gmail_user
        self.gmail_pwd = gmail_pwd

    def send(self, to, subject, text = "", attach = None):
        self.to = to
        self.subject = subject
        self.text = text
        self.attach = attach
        msg = MIMEMultipart()
        #################
        msg['From'] = self.gmail_user
        msg['To'] = self.to
        msg['Subject'] = self.subject
        #############
        if text != '' :
           msg.attach(MIMEText(self.text))
        ###############
        if attach != None :
           part = MIMEBase('application', 'octet-stream')
           part.set_payload(open(self.attach, 'rb').read())
           Encoders.encode_base64(part)
           part.add_header('Content-Disposition',
                   'attachment; filename = "%s"' % op.basename(self.attach))
           msg.attach(part)
        ##############
        mailServer = smtplib.SMTP("smtp.gmail.com", 587)
        mailServer.ehlo()
        mailServer.starttls()
        mailServer.ehlo()
        mailServer.login(self.gmail_user, self.gmail_pwd)
        mailServer.sendmail(self.gmail_user, self.to, msg.as_string())
        mailServer.close()
        
class Test(unittest.TestCase):
    gm = GMAIL()
    gm.send(to = 'gmalert67@gmail.com', subject = 'test gmail', text = "hello", attach = None)

if __name__ == '__main__':
    unittest.main()
