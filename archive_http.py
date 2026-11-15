"""Bounded HTTP retries for transient archive and upstream gateway failures."""
import math
import time
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from urllib3.exceptions import MaxRetryError, ResponseError


class BoundedRetry(Retry):
    """One shared retry budget and bounded server-directed waits.

    urllib3's allowed_methods gate does not cover connection failures. Guard
    increment too so non-idempotent methods are not retried after any failure.
    """
    def increment(self, method=None, url=None, response=None, error=None, _pool=None, _stacktrace=None):
        if method is not None and not self._is_method_retryable(method):
            reason=error or ResponseError('Retries disabled for this HTTP method')
            raise MaxRetryError(_pool, url, reason) from error
        return super().increment(method, url, response, error, _pool, _stacktrace)

    def sleep_for_retry(self, response):
        delay=self.get_retry_after(response)
        if delay:
            time.sleep(min(float(delay),4.))
            return True
        return False


def retry_session(*, total=3, backoff=0.5, pool_connections=4, pool_maxsize=8):
    """Finite retries for idempotent GET/HEAD requests (at most 11 attempts)."""
    if not isinstance(total, int) or not 0 <= total <= 10:
        raise ValueError('total retries must be between 0 and 10')
    if not math.isfinite(backoff) or backoff<0:raise ValueError('backoff must be finite and nonnegative')
    retry = BoundedRetry(total=total, connect=total, read=total, status=total,
                  backoff_factor=backoff, backoff_max=4, backoff_jitter=.25,
                  status_forcelist=(408, 425, 429, 500, 502, 503, 504),
                  allowed_methods=frozenset({'GET', 'HEAD'}),
                  respect_retry_after_header=True, raise_on_status=True)
    adapter = HTTPAdapter(max_retries=retry, pool_connections=pool_connections,
                          pool_maxsize=pool_maxsize)
    session = requests.Session()
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    return session
