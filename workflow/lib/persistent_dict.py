"""Generic persistent, concurrent dictionary-like facility."""

"""
Copied into nextstrain/ncov workflow from pytools version 2021.2.9 unchanged (except for this message).
https://github.com/inducer/pytools/blob/a0d9fdc33f459802b54a79c00a617fbc410e9d4d/pytools/persistent_dict.py.
See https://snakemake.readthedocs.io/en/stable/project_info/faq.html#i-want-to-pass-variables-between-rules-is-that-possible
for rationale.
"""

__copyright__ = """
Copyright (C) 2011,2014 Andreas Kloeckner
Copyright (C) 2017 Matt Wala
"""

__license__ = """
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import logging
import hashlib
import collections.abc as abc

# Removing this in 2020-12 broke a shocking amount of stuff, such as
# https://github.com/OP2/PyOP2/pull/605
# Bring it back for now to mitigate the breakage, however this is going
# away in 2021 at the latest.
new_hash = hashlib.sha256

import os
import shutil
import sys
import errno

logger = logging.getLogger(__name__)

__doc__ = """
Persistent Hashing and Persistent Dictionaries
==============================================

This module contains functionality that allows hashing with keys that remain
valid across interpreter invocations, unlike Python's built-in hashes.

This module also provides a disk-backed dictionary that uses persistent hashing.

.. autoexception:: NoSuchEntryError
.. autoexception:: ReadOnlyEntryError

.. autoexception:: CollisionWarning

.. autoclass:: KeyBuilder
.. autoclass:: PersistentDict
.. autoclass:: WriteOncePersistentDict
"""


def _make_dir_recursively(dir_):
    try:
        os.makedirs(dir_)
    except OSError as ex:
        from errno import EEXIST
        if ex.errno != EEXIST:
            raise


def update_checksum(checksum, obj):
    if isinstance(obj, str):
        checksum.update(obj.encode("utf8"))
    else:
        checksum.update(obj)


# {{{ cleanup managers

class CleanupBase:
    pass


class CleanupManager(CleanupBase):
    def __init__(self):
        self.cleanups = []

    def register(self, c):
        self.cleanups.insert(0, c)

    def clean_up(self):
        for c in self.cleanups:
            c.clean_up()

    def error_clean_up(self):
        for c in self.cleanups:
            c.error_clean_up()


class LockManager(CleanupBase):
    def __init__(self, cleanup_m, lock_file, stacklevel=0):
        self.lock_file = lock_file

        attempts = 0
        while True:
            try:
                self.fd = os.open(self.lock_file,
                        os.O_CREAT | os.O_WRONLY | os.O_EXCL)
                break
            except OSError:
                pass

            from time import sleep
            sleep(1)

            attempts += 1

            if attempts > 10:
                from warnings import warn
                warn("could not obtain lock -- "
                        f"delete '{self.lock_file}' if necessary",
                        stacklevel=1 + stacklevel)
            if attempts > 3 * 60:
                raise RuntimeError("waited more than three minutes "
                        f"on the lock file '{self.lock_file}' "
                        "-- something is wrong")

        cleanup_m.register(self)

    def clean_up(self):
        os.close(self.fd)
        os.unlink(self.lock_file)

    def error_clean_up(self):
        pass


class ItemDirManager(CleanupBase):
    def __init__(self, cleanup_m, path, delete_on_error):
        from os.path import isdir

        self.existed = isdir(path)
        self.path = path
        self.delete_on_error = delete_on_error

        cleanup_m.register(self)

    def reset(self):
        try:
            shutil.rmtree(self.path)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise

    def mkdir(self):
        from os import makedirs
        makedirs(self.path, exist_ok=True)

    def clean_up(self):
        pass

    def error_clean_up(self):
        if self.delete_on_error:
            self.reset()

# }}}


# {{{ key generation

class KeyBuilder:
    """A (stateless) object that computes hashes of objects fed to it. Subclassing
    this class permits customizing the computation of hash keys.

    .. automethod:: __call__
    .. automethod:: rec
    .. staticmethod:: new_hash()

        Return a new hash instance following the protocol of the ones
        from :mod:`hashlib`. This will permit switching to different
        hash algorithms in the future. Subclasses are expected to use
        this to create new hashes. Not doing so is deprecated and
        may stop working as early as 2022.

        .. versionadded:: 2021.2
    """

    # this exists so that we can (conceivably) switch algorithms at some point
    # down the road
    new_hash = hashlib.sha256

    def rec(self, key_hash, key):
        """
        :arg key_hash: the hash object to be updated with the hash of *key*.
        :arg key: the (immutable) Python object to be hashed.
        :returns: the updated *key_hash*

        .. versionchanged:: 2021.2

            Now returns the updated *key_hash*.
        """

        digest = None

        try:
            digest = key._pytools_persistent_hash_digest  # noqa pylint:disable=protected-access
        except AttributeError:
            pass

        if digest is None:
            try:
                method = key.update_persistent_hash
            except AttributeError:
                pass
            else:
                inner_key_hash = self.new_hash()
                method(inner_key_hash, self)
                digest = inner_key_hash.digest()

        if digest is None:
            tname = type(key).__name__
            method = None
            try:
                method = getattr(self, "update_for_"+tname)
            except AttributeError:
                # Handling numpy >= 1.20, for which
                # type(np.dtype("float32")) -> "dtype[float32]"
                if tname.startswith("dtype[") and "numpy" in sys.modules:
                    import numpy as np
                    if isinstance(key, np.dtype):
                        method = self.update_for_specific_dtype

            if method is not None:
                inner_key_hash = self.new_hash()
                method(inner_key_hash, key)
                digest = inner_key_hash.digest()

        if digest is None:
            raise TypeError(
                    f"unsupported type for persistent hash keying: {type(key)}")

        if not isinstance(key, type):
            try:
                key._pytools_persistent_hash_digest = digest   # noqa pylint:disable=protected-access
            except AttributeError:
                pass
            except TypeError:
                pass

        key_hash.update(digest)
        return key_hash

    def __call__(self, key):
        key_hash = self.new_hash()
        self.rec(key_hash, key)
        return key_hash.hexdigest()

    # {{{ updaters

    @staticmethod
    def update_for_int(key_hash, key):
        sz = 8
        while True:
            try:
                key_hash.update(key.to_bytes(sz, byteorder="little", signed=True))
                return
            except OverflowError:
                sz *= 2

    @staticmethod
    def update_for_bool(key_hash, key):
        key_hash.update(str(key).encode("utf8"))

    @staticmethod
    def update_for_float(key_hash, key):
        key_hash.update(key.hex().encode("utf8"))

    @staticmethod
    def update_for_str(key_hash, key):
        key_hash.update(key.encode("utf8"))

    @staticmethod
    def update_for_bytes(key_hash, key):
        key_hash.update(key)

    def update_for_tuple(self, key_hash, key):
        for obj_i in key:
            self.rec(key_hash, obj_i)

    def update_for_frozenset(self, key_hash, key):
        from pytools import unordered_hash

        unordered_hash(
            key_hash,
            (self.rec(self.new_hash(), key_i).digest() for key_i in key))

    @staticmethod
    def update_for_NoneType(key_hash, key):  # noqa
        del key
        key_hash.update(b"<None>")

    @staticmethod
    def update_for_dtype(key_hash, key):
        key_hash.update(key.str.encode("utf8"))

    # Handling numpy >= 1.20, for which
    # type(np.dtype("float32")) -> "dtype[float32]"
    # Introducing this method allows subclasses to specially handle all those
    # dtypes.
    @staticmethod
    def update_for_specific_dtype(key_hash, key):
        key_hash.update(key.str.encode("utf8"))

    # }}}

# }}}


# {{{ lru cache

class _LinkedList:
    """The list operates on nodes of the form [value, leftptr, rightpr]. To create a
    node of this form you can use `LinkedList.new_node().`

    Supports inserting at the left and deleting from an arbitrary location.
    """
    def __init__(self):
        self.count = 0
        self.head = None
        self.end = None

    @staticmethod
    def new_node(element):
        return [element, None, None]

    def __len__(self):
        return self.count

    def appendleft_node(self, node):
        self.count += 1

        if self.head is None:
            self.head = self.end = node
            return

        self.head[1] = node
        node[2] = self.head

        self.head = node

    def pop_node(self):
        end = self.end
        self.remove_node(end)
        return end

    def remove_node(self, node):
        self.count -= 1

        if self.head is self.end:
            assert node is self.head
            self.head = self.end = None
            return

        left = node[1]
        right = node[2]

        if left is None:
            self.head = right
        else:
            left[2] = right

        if right is None:
            self.end = left
        else:
            right[1] = left

        node[1] = node[2] = None


class _LRUCache(abc.MutableMapping):
    """A mapping that keeps at most *maxsize* items with an LRU replacement policy.
    """
    def __init__(self, maxsize):
        self.lru_order = _LinkedList()
        self.maxsize = maxsize
        self.cache = {}

    def __delitem__(self, item):
        node = self.cache[item]
        self.lru_order.remove_node(node)
        del self.cache[item]

    def __getitem__(self, item):
        node = self.cache[item]
        self.lru_order.remove_node(node)
        self.lru_order.appendleft_node(node)
        # A linked list node contains a tuple of the form (item, value).
        return node[0][1]

    def __contains__(self, item):
        return item in self.cache

    def __iter__(self):
        return iter(self.cache)

    def __len__(self):
        return len(self.cache)

    def clear(self):
        self.cache.clear()
        self.lru_order = _LinkedList()

    def __setitem__(self, item, value):
        if self.maxsize < 1:
            return

        try:
            node = self.cache[item]
            self.lru_order.remove_node(node)
        except KeyError:
            if len(self.lru_order) >= self.maxsize:
                # Make room for new elements.
                end_node = self.lru_order.pop_node()
                del self.cache[end_node[0][0]]

            node = self.lru_order.new_node((item, value))
            self.cache[item] = node

        self.lru_order.appendleft_node(node)

        assert len(self.cache) == len(self.lru_order), \
                (len(self.cache), len(self.lru_order))
        assert len(self.lru_order) <= self.maxsize

# }}}


# {{{ top-level

class NoSuchEntryError(KeyError):
    pass


class ReadOnlyEntryError(KeyError):
    pass


class CollisionWarning(UserWarning):
    pass


class _PersistentDictBase:
    def __init__(self, identifier, key_builder=None, container_dir=None):
        self.identifier = identifier

        if key_builder is None:
            key_builder = KeyBuilder()

        self.key_builder = key_builder

        from os.path import join
        if container_dir is None:
            try:
                import platformdirs as appdirs
            except ImportError:
                import appdirs

            container_dir = join(
                    appdirs.user_cache_dir("pytools", "pytools"),
                    "pdict-v4-{}-py{}".format(
                        identifier,
                        ".".join(str(i) for i in sys.version_info)))

        self.container_dir = container_dir

        self._make_container_dir()

    @staticmethod
    def _warn(msg, category=UserWarning, stacklevel=0):
        from warnings import warn
        warn(msg, category, stacklevel=1 + stacklevel)

    def store_if_not_present(self, key, value, _stacklevel=0):
        self.store(key, value, _skip_if_present=True, _stacklevel=1 + _stacklevel)

    def store(self, key, value, _skip_if_present=False, _stacklevel=0):
        raise NotImplementedError()

    def fetch(self, key, _stacklevel=0):
        raise NotImplementedError()

    @staticmethod
    def _read(path):
        from pickle import load
        with open(path, "rb") as inf:
            return load(inf)

    @staticmethod
    def _write(path, value):
        from pickle import dump, HIGHEST_PROTOCOL
        with open(path, "wb") as outf:
            dump(value, outf, protocol=HIGHEST_PROTOCOL)

    def _item_dir(self, hexdigest_key):
        from os.path import join
        # Some file systems limit the number of directories in a directory.
        # For ext4, that limit appears to be 64K for example.
        # This doesn't solve that problem, but it makes it much less likely

        return join(self.container_dir,
                hexdigest_key[:3],
                hexdigest_key[3:6],
                hexdigest_key[6:])

    def _key_file(self, hexdigest_key):
        from os.path import join
        return join(self._item_dir(hexdigest_key), "key")

    def _contents_file(self, hexdigest_key):
        from os.path import join
        return join(self._item_dir(hexdigest_key), "contents")

    def _lock_file(self, hexdigest_key):
        from os.path import join
        return join(self.container_dir, str(hexdigest_key) + ".lock")

    def _make_container_dir(self):
        _make_dir_recursively(self.container_dir)

    def _collision_check(self, key, stored_key, _stacklevel):
        if stored_key != key:
            # Key collision, oh well.
            self._warn(f"{self.identifier}: key collision in cache at "
                    f"'{self.container_dir}' -- these are sufficiently unlikely "
                    "that they're often indicative of a broken hash key "
                    "implementation (that is not considering some elements "
                    "relevant for equality comparison)",
                    CollisionWarning,
                    1 + _stacklevel)

            # This is here so we can step through equality comparison to
            # see what is actually non-equal.
            stored_key == key  # pylint:disable=pointless-statement  # noqa: B015
            raise NoSuchEntryError(key)

    def __getitem__(self, key):
        return self.fetch(key, _stacklevel=1)

    def __setitem__(self, key, value):
        self.store(key, value, _stacklevel=1)

    def clear(self):
        try:
            shutil.rmtree(self.container_dir)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise

        self._make_container_dir()


class WriteOncePersistentDict(_PersistentDictBase):
    """A concurrent disk-backed dictionary that disallows overwriting/deletion.

    Compared with :class:`PersistentDict`, this class has faster
    retrieval times.

    .. automethod:: __init__
    .. automethod:: __getitem__
    .. automethod:: __setitem__
    .. automethod:: clear
    .. automethod:: store
    .. automethod:: store_if_not_present
    .. automethod:: fetch
    """
    def __init__(self, identifier, key_builder=None, container_dir=None,
             in_mem_cache_size=256):
        """
        :arg identifier: a file-name-compatible string identifying this
            dictionary
        :arg key_builder: a subclass of :class:`KeyBuilder`
        :arg in_mem_cache_size: retain an in-memory cache of up to
            *in_mem_cache_size* items
        """
        _PersistentDictBase.__init__(self, identifier, key_builder, container_dir)
        self._cache = _LRUCache(in_mem_cache_size)

    def _spin_until_removed(self, lock_file, stacklevel):
        from os.path import exists

        attempts = 0
        while exists(lock_file):
            from time import sleep
            sleep(1)

            attempts += 1

            if attempts > 10:
                self._warn(
                        f"waiting until unlocked--delete '{lock_file}' if necessary",
                        stacklevel=1 + stacklevel)

            if attempts > 3 * 60:
                raise RuntimeError("waited more than three minutes "
                        f"on the lock file '{lock_file}'"
                        "--something is wrong")

    def store(self, key, value, _skip_if_present=False, _stacklevel=0):
        hexdigest_key = self.key_builder(key)

        cleanup_m = CleanupManager()
        try:
            try:
                LockManager(cleanup_m, self._lock_file(hexdigest_key),
                        1 + _stacklevel)
                item_dir_m = ItemDirManager(
                        cleanup_m, self._item_dir(hexdigest_key),
                        delete_on_error=False)

                if item_dir_m.existed:
                    if _skip_if_present:
                        return
                    raise ReadOnlyEntryError(key)

                item_dir_m.mkdir()

                key_path = self._key_file(hexdigest_key)
                value_path = self._contents_file(hexdigest_key)

                self._write(value_path, value)
                self._write(key_path, key)

                logger.debug("%s: disk cache store [key=%s]",
                        self.identifier, hexdigest_key)
            except Exception:
                cleanup_m.error_clean_up()
                raise
        finally:
            cleanup_m.clean_up()

    def fetch(self, key, _stacklevel=0):
        hexdigest_key = self.key_builder(key)

        # {{{ in memory cache

        try:
            stored_key, stored_value = self._cache[hexdigest_key]
        except KeyError:
            pass
        else:
            logger.debug("%s: in mem cache hit [key=%s]",
                    self.identifier, hexdigest_key)
            self._collision_check(key, stored_key, 1 + _stacklevel)
            return stored_value

        # }}}

        # {{{ check path exists and is unlocked

        item_dir = self._item_dir(hexdigest_key)

        from os.path import isdir
        if not isdir(item_dir):
            logger.debug("%s: disk cache miss [key=%s]",
                    self.identifier, hexdigest_key)
            raise NoSuchEntryError(key)

        lock_file = self._lock_file(hexdigest_key)
        self._spin_until_removed(lock_file, 1 + _stacklevel)

        # }}}

        key_file = self._key_file(hexdigest_key)
        contents_file = self._contents_file(hexdigest_key)

        # Note: Unlike PersistentDict, this doesn't autodelete invalid entires,
        # because that would lead to a race condition.

        # {{{ load key file and do equality check

        try:
            read_key = self._read(key_file)
        except Exception as e:
            self._warn(f"{type(self).__name__}({self.identifier}) "
                    f"encountered an invalid key file for key {hexdigest_key}. "
                    f"Remove the directory '{item_dir}' if necessary. "
                    f"(caught: {type(e).__name__}: {e})",
                    stacklevel=1 + _stacklevel)
            raise NoSuchEntryError(key)

        self._collision_check(key, read_key, 1 + _stacklevel)

        # }}}

        logger.debug("%s: disk cache hit [key=%s]",
                self.identifier, hexdigest_key)

        # {{{ load contents

        try:
            read_contents = self._read(contents_file)
        except Exception:
            self._warn(f"{type(self).__name__}({self.identifier}) "
                    f"encountered an invalid key file for key {hexdigest_key}. "
                    f"Remove the directory '{item_dir}' if necessary.",
                    stacklevel=1 + _stacklevel)
            raise NoSuchEntryError(key)

        # }}}

        self._cache[hexdigest_key] = (key, read_contents)
        return read_contents

    def clear(self):
        _PersistentDictBase.clear(self)
        self._cache.clear()


class PersistentDict(_PersistentDictBase):
    """A concurrent disk-backed dictionary.

    .. automethod:: __init__
    .. automethod:: __getitem__
    .. automethod:: __setitem__
    .. automethod:: __delitem__
    .. automethod:: clear
    .. automethod:: store
    .. automethod:: store_if_not_present
    .. automethod:: fetch
    .. automethod:: remove
    """
    def __init__(self, identifier, key_builder=None, container_dir=None):
        """
        :arg identifier: a file-name-compatible string identifying this
            dictionary
        :arg key_builder: a subclass of :class:`KeyBuilder`
        """
        _PersistentDictBase.__init__(self, identifier, key_builder, container_dir)

    def store(self, key, value, _skip_if_present=False, _stacklevel=0):
        hexdigest_key = self.key_builder(key)

        cleanup_m = CleanupManager()
        try:
            try:
                LockManager(cleanup_m, self._lock_file(hexdigest_key),
                        1 + _stacklevel)
                item_dir_m = ItemDirManager(
                        cleanup_m, self._item_dir(hexdigest_key),
                        delete_on_error=True)

                if item_dir_m.existed:
                    if _skip_if_present:
                        return
                    item_dir_m.reset()

                item_dir_m.mkdir()

                key_path = self._key_file(hexdigest_key)
                value_path = self._contents_file(hexdigest_key)

                self._write(value_path, value)
                self._write(key_path, key)

                logger.debug("%s: cache store [key=%s]",
                        self.identifier, hexdigest_key)
            except Exception:
                cleanup_m.error_clean_up()
                raise
        finally:
            cleanup_m.clean_up()

    def fetch(self, key, _stacklevel=0):
        hexdigest_key = self.key_builder(key)
        item_dir = self._item_dir(hexdigest_key)

        from os.path import isdir
        if not isdir(item_dir):
            logger.debug("%s: cache miss [key=%s]",
                    self.identifier, hexdigest_key)
            raise NoSuchEntryError(key)

        cleanup_m = CleanupManager()
        try:
            try:
                LockManager(cleanup_m, self._lock_file(hexdigest_key),
                        1 + _stacklevel)
                item_dir_m = ItemDirManager(
                        cleanup_m, item_dir, delete_on_error=False)

                key_path = self._key_file(hexdigest_key)
                value_path = self._contents_file(hexdigest_key)

                # {{{ load key

                try:
                    read_key = self._read(key_path)
                except Exception:
                    item_dir_m.reset()
                    self._warn(f"{type(self).__name__}({self.identifier}) "
                            "encountered an invalid key file for key "
                            f"{hexdigest_key}. Entry deleted.",
                            stacklevel=1 + _stacklevel)
                    raise NoSuchEntryError(key)

                self._collision_check(key, read_key, 1 + _stacklevel)

                # }}}

                logger.debug("%s: cache hit [key=%s]",
                        self.identifier, hexdigest_key)

                # {{{ load value

                try:
                    read_contents = self._read(value_path)
                except Exception:
                    item_dir_m.reset()
                    self._warn(f"{type(self).__name__}({self.identifier}) "
                            "encountered an invalid key file for key "
                            f"{hexdigest_key}. Entry deleted.",
                            stacklevel=1 + _stacklevel)
                    raise NoSuchEntryError(key)

                return read_contents

                # }}}

            except Exception:
                cleanup_m.error_clean_up()
                raise
        finally:
            cleanup_m.clean_up()

    def remove(self, key, _stacklevel=0):
        hexdigest_key = self.key_builder(key)

        item_dir = self._item_dir(hexdigest_key)
        from os.path import isdir
        if not isdir(item_dir):
            raise NoSuchEntryError(key)

        cleanup_m = CleanupManager()
        try:
            try:
                LockManager(cleanup_m, self._lock_file(hexdigest_key),
                        1 + _stacklevel)
                item_dir_m = ItemDirManager(
                        cleanup_m, item_dir, delete_on_error=False)
                key_file = self._key_file(hexdigest_key)

                # {{{ load key

                try:
                    read_key = self._read(key_file)
                except Exception:
                    item_dir_m.reset()
                    self._warn(f"{type(self).__name__}({self.identifier}) "
                            "encountered an invalid key file for key "
                            f"{hexdigest_key}. Entry deleted",
                            stacklevel=1 + _stacklevel)
                    raise NoSuchEntryError(key)

                self._collision_check(key, read_key, 1 + _stacklevel)

                # }}}

                item_dir_m.reset()

            except Exception:
                cleanup_m.error_clean_up()
                raise
        finally:
            cleanup_m.clean_up()

    def __delitem__(self, key):
        self.remove(key, _stacklevel=1)

# }}}

# vim: foldmethod=marker
