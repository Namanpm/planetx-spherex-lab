"""Transport-level retry checks use a real loopback server, never NASA."""
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from threading import Thread
import pytest
import requests
from urllib3.exceptions import MaxRetryError, NewConnectionError
from urllib3.response import HTTPResponse
from archive_http import BoundedRetry, retry_session


@contextmanager
def server(statuses):
    calls=[]
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            calls.append(self.command)
            code=statuses[min(len(calls)-1,len(statuses)-1)]
            self.send_response(code);self.send_header('Content-Length','2');self.end_headers();self.wfile.write(b'ok')
        do_POST=do_GET
        def log_message(self,*args):pass
    http=ThreadingHTTPServer(('127.0.0.1',0),Handler)
    thread=Thread(target=http.serve_forever,kwargs={'poll_interval':.01},daemon=True);thread.start()
    try:yield f'http://127.0.0.1:{http.server_address[1]}/query',calls
    finally:http.shutdown();http.server_close();thread.join(timeout=2)


@pytest.fixture
def no_delay(monkeypatch):
    sleeps=[]
    monkeypatch.setattr(BoundedRetry,'sleep',lambda self,response=None:sleeps.append(self.get_backoff_time()))
    return sleeps


def test_real_http_502_503_recovery_uses_combined_four_attempt_budget(no_delay):
    with server([502,503,502,200]) as (url,calls),retry_session() as client:
        response=client.get(url,timeout=(1,1))
        assert response.status_code==200
        assert calls==['GET']*4
        assert len(no_delay)==3
        assert all(0<=delay<=4 for delay in no_delay)


def test_real_http_retry_exhaustion_stops_after_four_attempts(no_delay):
    with server([502,503,502,503,200]) as (url,calls),retry_session() as client:
        with pytest.raises(requests.exceptions.RetryError):client.get(url,timeout=(1,1))
        assert len(calls)==4


def test_permanent_400_is_not_retried(no_delay):
    with server([400,200]) as (url,calls),retry_session() as client:
        response=client.get(url,timeout=(1,1))
        with pytest.raises(requests.HTTPError):response.raise_for_status()
        assert calls==['GET']


def test_post_is_never_retried_after_gateway_failure(no_delay):
    with server([503,200]) as (url,calls),retry_session() as client:
        assert client.post(url,data='request',timeout=(1,1)).status_code==503
        assert calls==['POST']


def test_post_connection_failure_cannot_consume_retry_budget():
    retry=retry_session().get_adapter('https://').max_retries
    with pytest.raises(MaxRetryError):
        retry.increment('POST','https://example.invalid',error=NewConnectionError(None,'simulated connection failure'))


def test_retry_after_is_clamped_to_four_seconds(monkeypatch):
    sleeps=[]
    monkeypatch.setattr('archive_http.time.sleep',sleeps.append)
    retry=retry_session().get_adapter('https://').max_retries
    assert retry.sleep_for_retry(HTTPResponse(headers={'Retry-After':'9999999'}))
    assert sleeps==[4.]


def test_exponential_jitter_is_bounded_and_survives_policy_cloning():
    retry=retry_session(backoff=.5).get_adapter('https://').max_retries
    retry=retry.increment('GET','https://example.invalid',response=HTTPResponse(status=502))
    retry=retry.increment('GET','https://example.invalid',response=HTTPResponse(status=503))
    assert retry.total==1
    assert retry.backoff_jitter==.25
    assert 1.<=retry.get_backoff_time()<=1.25
    high=retry_session(backoff=100.).get_adapter('https://').max_retries
    for _ in range(2):high=high.increment('GET','https://example.invalid',response=HTTPResponse(status=502))
    assert high.get_backoff_time()==4.


@pytest.mark.parametrize('total,backoff',[(11,.5),(-1,.5),(1,-1),(1,float('inf'))])
def test_invalid_retry_limits_rejected(total,backoff):
    with pytest.raises(ValueError):retry_session(total=total,backoff=backoff)
